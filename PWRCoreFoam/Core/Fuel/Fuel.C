/*----------------------------------------------------------------------------*\
| Class:                                                                       |
|     Foam::Fuel                                                               |
\*----------------------------------------------------------------------------*/


#include "Fuel.H"

#include "fvCFD.H"
#include "rhoThermo.H"
#include "rodFieldsFwd.H"
#include "subChannelMesh.H"


//****************************************************************************//

namespace Foam
{

defineTypeNameAndDebug(Fuel, 0);

}


/*----------------------------------------------------------------------------*\
|                       Constructors                                           |
\*----------------------------------------------------------------------------*/

Foam::Fuel::Fuel
(
  const subChannelMesh & mesh,
  const volVectorField & U,
  const rhoThermo & thermo
)
:
  fuelDict_
  (
    IOobject
    (
      "ReactorCore",
      mesh.time().constant(),
      mesh,
      IOobject::MUST_READ_IF_MODIFIED,
      IOobject::NO_WRITE
    )
  ),

  mesh_(mesh),

  U_(U),

  thermo_(thermo),

  Cp_
  (
    IOobject
    (
      "Cp_",
      mesh.time().timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    thermo.Cp()()
  ),

  kappa_
  (
    IOobject
    (
      "kappa_",
      mesh.time().timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    thermo.kappa()()
  ),

  rodLinePowerDensity_
  (
    IOobject
    (
      "rodLinearPowerDensity",
      mesh.time().timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rodLinearPowerDensity", dimensionSet(1, 1, -3, 0, 0, 0, 0), Foam::scalar(0.0))
  ),

  pelletCenterTemperature_
  (
    IOobject
    (
      "pelletCenterTemperature",
      mesh.time().timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("pelletCenterTemperature", dimensionSet(0, 0, 0, 1, 0, 0, 0), Foam::scalar(0.0))
  ),

  cladSurfaceTemperature_
  (
    IOobject
    (
      "cladSurfaceTemperature",
      mesh.time().timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("cladSurfaceTemperature", dimensionSet(0, 0, 0, 1, 0, 0, 0), Foam::scalar(0.0))
  ),

  powerTypeTable_(0),

  Re_for_rod_
  (
    IOobject
    (
      "Re_for_rod",
      mesh.time().timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Re_for_rod_", dimensionSet(0, 0, 0, 1, 0, 0, 0), Foam::scalar(0.0))
  ),

  h_for_rod_
  (
    IOobject
    (
      "h_for_rod",
      mesh.time().timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("h_for_rod_", dimensionSet(0, 0, 0, 1, 0, 0, 0), Foam::scalar(0.0))
  ),

  U_time_rho_
  (
    IOobject
    (
      "U_time_rho",
      mesh.time().timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("U_time_rho", dimensionSet(0, 0, 0, 1, 0, 0, 0), Foam::scalar(0.0))
  ),

  T_coolant_
  (
    IOobject
    (
      "T_coolant",
      mesh.time().timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("T_coolant", dimensionSet(0, 0, 0, 1, 0, 0, 0), Foam::scalar(0.0))
  )

{
  set_fuel_power();
}


/*----------------------------------------------------------------------------*\
|                       Destructor                                             |
\*----------------------------------------------------------------------------*/
Foam::Fuel::~Fuel()
{}


/*----------------------------------------------------------------------------*\
|                       Member functions                                       |
\*----------------------------------------------------------------------------*/

// Private member functions

void
Foam::Fuel::set_uniform_power_distribution(scalar uniform_line_power_density)
{
  const labelList & rodEdgesList_ = mesh_.get_rodEdgesList();

  if (!(rodLinePowerDensity_.size() == rodEdgesList_.size()))
  {
    FatalErrorIn("void Foam::Fuel::set_uniform_power_distribution()")
              << "rodLinePowerDensity_ size does not equal to rodEdgesList_ size"
              << exit(FatalError);
  }

  forAll(rodEdgesList_, i)
  {
    rodLinePowerDensity_[i] = uniform_line_power_density;
  }
}


void
Foam::Fuel::set_cos_power_distribution(scalar average_line_power_density)
{
  const labelList & rodEdgesList_ = mesh_.get_rodEdgesList();

  if (!(rodLinePowerDensity_.size() == rodEdgesList_.size()))
  {
    FatalErrorIn("void Foam::Fuel::set_cos_power_distribution()")
              << "rodLinePowerDensity_ size does not equal to rodEdgesList_ size"
              << exit(FatalError);
  }

  scalar pi = 3.14159265;
  scalar power_max = average_line_power_density * pi / 2.0;

  scalar z_max = rodLinePowerDensity_.get_z_max();
  scalar z_max_global = returnReduce(z_max, maxOp<scalar>(), Pstream::msgType(), UPstream::worldComm);

  scalar z_min = rodLinePowerDensity_.get_z_min();
  scalar z_min_global = returnReduce(z_min, minOp<scalar>(), Pstream::msgType(), UPstream::worldComm);

  scalar core_height = z_max_global - z_min_global;

  forAll(rodEdgesList_, i)
  {
    label edge_i = rodEdgesList_[i];

    scalar local_power = 0.0;

    scalar z_1 = mesh_.get_rod_edge_z_min(edge_i);
    scalar z_2 = mesh_.get_rod_edge_z_max(edge_i);

    local_power = -1.0 * power_max * core_height / pi *
                  (cos(z_2 * pi / core_height) - cos(z_1 * pi / core_height)) /
                  (z_2 - z_1);

    rodLinePowerDensity_[i] = local_power;
  }
}


void
Foam::Fuel::set_polynominal_power_distribution(scalar average_line_power_density)
{
  const labelList & rodEdgesList_ = mesh_.get_rodEdgesList();

  if (!(rodLinePowerDensity_.size() == rodEdgesList_.size()))
  {
    FatalErrorIn("void Foam::Fuel::set_polynominal_power_distribution()")
              << "rodLinePowerDensity_ size does not equal to rodEdgesList_ size"
              << exit(FatalError);
  }

  scalarList power_coeffs(fuelDict_.lookup("PowerCoeffs"));

  forAll(rodEdgesList_, i)
  {
    label edge_i = rodEdgesList_[i];

    scalar local_power = average_line_power_density;

    scalar z_1 = mesh_.get_rod_edge_z_min(edge_i);
    scalar z_2 = mesh_.get_rod_edge_z_max(edge_i);

    local_power *= calc_integral_of_polynominal(z_1, z_2, power_coeffs) / (z_2 - z_1);

    rodLinePowerDensity_[i] = local_power;
  }
}


scalar
Foam::Fuel::calc_value_of_polynominal
(
  const scalar & x,
  const scalarList & polynominal_coeffs
)
{
  scalar retval = 0.0;

  forAll(polynominal_coeffs, i)
  {
    retval += polynominal_coeffs[i] *
              ::pow(x, scalar(i));
  }

  return retval;
}


scalar
Foam::Fuel::calc_integral_of_polynominal
(
  const scalar & x_1,
  const scalar & x_2,
  const scalarList & polynominal_coeffs
)
{
  scalar integral_1 = 0.0;
  scalar integral_2 = 0.0;

  forAll(polynominal_coeffs, i)
  {
    integral_1 += polynominal_coeffs[i] / scalar(i + 1.0) *
                  ::pow(x_1, scalar(i+1.0));

    integral_2 += polynominal_coeffs[i] / scalar(i + 1.0) *
                  ::pow(x_2, scalar(i+1.0));
  }

  return (integral_2 - integral_1);
}


void
Foam::Fuel::set_horizontal_power_distribution()
{
  forAll(rodLinePowerDensity_, i)
  {
    label rod_x = mesh_.get_rod_rod_index_x_y()[i][0];
    label rod_y = mesh_.get_rod_rod_index_x_y()[i][1];

    label assembly_x = mesh_.get_rod_assembly_index_x_y()[i][0];
    label assembly_y = mesh_.get_rod_assembly_index_x_y()[i][1];

    scalar factor = mesh_.get_rod_pin_power_distribution_factor()[rod_y-1][rod_x-1];
    scalar factor_assembly = mesh_.get_core_assembly_power_distribution_factor()[assembly_y-1][assembly_x-1];

    rodLinePowerDensity_[i] *= factor*factor_assembly;
  }
}


// Public member functions

void
Foam::Fuel::set_fuel_power()
{
  word power_type_(fuelDict_.lookup("FuelAxialPowerType"));
  label power_(readScalar(fuelDict_.lookup("FuelAverageLinePowerDensity")));

  if (powerTypeTable_.empty())
  {
    powerTypeTable_.insert("uniform",     1);
    powerTypeTable_.insert("cosin",       2);
    powerTypeTable_.insert("polynominal", 3);
  }

  switch(powerTypeTable_[power_type_])
  {
    case 1:
      set_uniform_power_distribution(power_);
      break;
    case 2:
      set_cos_power_distribution(power_);
      break;
    case 3:
      set_polynominal_power_distribution(power_);
      break;
    default:
      FatalErrorIn("void Foam::Fuel::set_fuel_power()")
                << "The FuelAxialPowerType option is not correct !"
                << exit(FatalError);
  }

  set_horizontal_power_distribution();
}


void
Foam::Fuel::update_fuel_power()
{

}


void
Foam::Fuel::update_clad_temperature()
{
  scalar pi = 3.14159;
  scalar rod_surface_area = 0.0;
  scalar diameter_o = readScalar(fuelDict_.lookup("RodOuterDiameter"));
  scalar pitch = readScalar(fuelDict_.lookup("pitch"));

  //scalar D_e = 4.0 * (pitch * pitch - pi * ::pow(diameter_o, 2.0) / 4.0) / (pi * diameter_o);
  scalar D_e = diameter_o;

  Cp_ = thermo_.Cp()();
  kappa_ = thermo_.kappa()();

  forAll(cladSurfaceTemperature_, i)
  {
    label edge_i = mesh_.get_rodEdgesList()[i];

    scalar edge_length = mesh_.get_rodEdgesLengthList()[i];
    rod_surface_area = pi * D_e * edge_length;
    scalar q_flux = rodLinePowerDensity_[i] * edge_length / rod_surface_area;

    scalar u_z = get_avg_U_z_around_rod(edge_i);
    scalar T_i = get_avg_T_around_rod(edge_i);
    scalar rho_i = get_avg_rho_around_rod(edge_i);
    scalar mu_i = get_avg_mu_around_rod(edge_i);
    scalar Cp_i = get_avg_Cp_around_rod(edge_i);
    scalar kappa_i = get_avg_kappa_around_rod(edge_i);

    scalar Pr_i = Cp_i * mu_i / kappa_i;
    scalar Re_i = u_z * D_e * rho_i / mu_i;
    scalar Nu_i = 0.023 * ::pow(Re_i, 0.8) * :: pow(Pr_i, 0.4);

    scalar h_i = Nu_i * kappa_i / D_e;
    h_i = 0.667 * h_i;

    cladSurfaceTemperature_[i] = T_i + q_flux / h_i;

    Re_for_rod_[i] = Re_i;
    h_for_rod_[i] = h_i;
    U_time_rho_[i] = u_z * rho_i;
    T_coolant_[i] = T_i;
  }
}


void
Foam::Fuel::update_pellet_temperature()
{
  scalar pi = 3.14159;
  forAll(pelletCenterTemperature_, i)
  {
    label edge_i = mesh_.get_rodEdgesList()[i];

    scalar edge_length = mesh_.get_rodEdgesLengthList()[i];

    scalar T_clad = cladSurfaceTemperature_[i];

    scalar q_l = rodLinePowerDensity_[i];

    scalar kappa_pellet = 4.5;

    //pelletCenterTemperature_[i] = T_clad + q_l / 4.0 / pi / kappa_pellet;
    pelletCenterTemperature_[i] = T_clad + q_l / 4.8 / pi / kappa_pellet;
  }
}


void
Foam::Fuel::update()
{
  update_fuel_power();

  update_clad_temperature();

  update_pellet_temperature();
}


void
Foam::Fuel::write()
{
  rodLinePowerDensity_.writeRodField();

  pelletCenterTemperature_.writeRodField();

  cladSurfaceTemperature_.writeRodField();

  Re_for_rod_.writeRodField();
  h_for_rod_.writeRodField();
  U_time_rho_.writeRodField();
  T_coolant_.writeRodField();
}


