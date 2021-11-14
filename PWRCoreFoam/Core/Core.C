/*----------------------------------------------------------------------------*\
| Class:                                                                       |
|     Foam::Core                                                               |
\*----------------------------------------------------------------------------*/


#include "Core.H"

#include "Fuel.H"
#include "fvCFD.H"
#include "IOdictionary.H"
#include "rhoThermo.H"
#include "subChannelMesh.H"


//****************************************************************************//

namespace Foam
{

defineTypeNameAndDebug(Core, 0);

}


/*----------------------------------------------------------------------------*\
|                       Constructors                                           |
\*----------------------------------------------------------------------------*/

Foam::Core::Core
(
  const subChannelMesh & mesh,
  const volVectorField & U,
  const rhoThermo & thermo
)
:
  Fuel(mesh, U, thermo),

  coreDict_
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

  momentumSource_
  (
    IOobject
    (
      "momentumSource",
      mesh.time().timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("momentumSource", dimensionSet(1, -2, -2, 0, 0, 0, 0), Foam::vector(0.0, 0.0, 0.0))
  ),

  energySource_
  (
    IOobject
    (
      "energySource",
      mesh.time().timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("energySource", dimensionSet(1, -1, -3, 0, 0, 0, 0), Foam::scalar(0.0))
  ),

  Re_
  (
    IOobject
    (
      "Re",
      mesh.time().timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Re", dimensionSet(0, 0, 0, 0, 0, 0, 0), Foam::scalar(0.0))
  ),

  phi_mix_
  (
    IOobject
    (
      "phi_mix",
      mesh.time().timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("phi_mix", dimensionSet(1, 0, -1, 0, 0), Foam::scalar(0.0))
  ),

  beta_(mesh.mixing_coefficient()),

  test_variable_1(&energySource_)
{
  update();

  write();
}


/*----------------------------------------------------------------------------*\
|                       Destructor                                             |
\*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*\
|                       Member functions                                       |
\*----------------------------------------------------------------------------*/

void
Foam::Core::calc_Re()
{
  forAll(Re_, cell_i)
  {
    scalar u_z = mag(U_[cell_i].z());
    scalar D_e = mesh_.get_hydraulic_diameter()[cell_i];
    scalar rho_i = thermo_.rho()()[cell_i];
    scalar mu_i = thermo_.mu()[cell_i];

    Re_[cell_i] = rho_i * u_z * D_e / mu_i;
  }
}


void
Foam::Core::update_momentum_source()
{
  forAll(momentumSource_, cell_i)
  {
    scalar u_z = U_[cell_i].z();
    scalar D_e = mesh_.get_hydraulic_diameter()[cell_i];
    scalar rho_i = thermo_.rho()()[cell_i];

    scalar Re_i = Re_[cell_i];
    scalar f_i = 0.0;

    if (Re_i <= 2300.0)
    {
      f_i = 64.0 / Re_i;
    }
    else
    {
      f_i = 0.0055 + 0.55 * ::pow(Re_i, -1.0/3.0);
    }

    momentumSource_[cell_i].z() = -1.0 / 2.0 * f_i * rho_i * mag(u_z) * u_z / D_e;

    label core_type_i = mesh_.get_core_type()[cell_i];

    if (core_type_i == 2)
    {
      momentumSource_[cell_i].z() += -1.0 / 2.0 * 1.04 * rho_i * mag(u_z) * u_z;
    }
    else if (core_type_i == 3)
    {
      momentumSource_[cell_i].z() += -1.0 / 2.0 * 1.61 * rho_i * mag(u_z) * u_z;
    }
  }
}


void
Foam::Core::calc_phi_mix()
{
  labelList gap_faces = mesh_.get_gapList();
  labelList owner_cells = mesh_.owner();
  labelList neighbour_cells = mesh_.neighbour();

  surfaceScalarField * magSfPtr_ = mesh_.get_magSfPtr();
  surfaceScalarField & faces_area_mag = *magSfPtr_;

  volScalarField rho = thermo_.rho()();

  forAll(gap_faces, gap_i)
  {
    label face_i = gap_faces[gap_i];
    label owner_i = owner_cells[face_i];
    label neighbour_i = neighbour_cells[face_i];

    scalar G_avg = 0.5 * (rho[owner_i]*U_[owner_i][2] + rho[neighbour_i]*U_[neighbour_i][2]);

    scalar W_mix = beta_ * G_avg * faces_area_mag[face_i];

    phi_mix_[face_i] = W_mix;
  }
}


void
Foam::Core::update_energy_source()
{
  forAll(energySource_, i)
  {
    energySource_[i] = 0.0;
  }

  forAll(get_rod_linear_power_density(), i)
  {
    label edge_i = mesh_.get_rodEdgesList()[i];
    scalar edge_length = mesh_.get_rodEdgesLengthList()[i];

    scalar rod_edge_total_power = edge_length * get_rod_linear_power_density()[i];

    forAll(mesh_.edgeCells()[edge_i], j)
    {
      label cell_j = mesh_.edgeCells()[edge_i][j];
      scalar cell_volume = mesh_.V()[cell_j];
      energySource_[cell_j] += rod_edge_total_power / 4.0 / cell_volume;
    }
  }
}


volScalarField *
Foam::Core::test_func_1() const
{
  return test_variable_1;
}


void
Foam::Core::update()
{
  Fuel::update();

  calc_Re();

  update_momentum_source();

  calc_phi_mix();

  update_energy_source();
}


void
Foam::Core::write()
{
  momentumSource_.write();

  energySource_.write();

  Re_.write();

  phi_mix_.write();

  Fuel::write();
}


