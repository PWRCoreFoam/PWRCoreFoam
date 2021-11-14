/*----------------------------------------------------------------------------*\
| Class:                                                                       |
|     Foam::Fuel                                                               |
\*----------------------------------------------------------------------------*/


//****************************************************************************//

namespace Foam
{

/*----------------------------------------------------------------------------*\
|                       Member functions                                       |
\*----------------------------------------------------------------------------*/

inline const Foam::IOdictionary & Foam::Fuel::get_fuel_dict()
{
  return fuelDict_;
}


inline const Foam::IOdictionary & Foam::Fuel::get_fuel_dict() const
{
  return fuelDict_;
}


inline Foam::rodScalarField & Foam::Fuel::get_rod_linear_power_density()
{
  return rodLinePowerDensity_;
}


inline const Foam::rodScalarField & Foam::Fuel::get_rod_linear_power_density() const
{
  return rodLinePowerDensity_;
}


inline Foam::rodScalarField & Foam::Fuel::get_pellet_center_temperature()
{
  return pelletCenterTemperature_;
}


inline const Foam::rodScalarField & Foam::Fuel::get_pellet_center_temperature() const
{
  return pelletCenterTemperature_;
}

inline Foam::rodScalarField & Foam::Fuel::get_clad_surface_temperature()
{
  return cladSurfaceTemperature_;
}


inline const Foam::rodScalarField & Foam::Fuel::get_clad_surface_temperature() const
{
  return cladSurfaceTemperature_;
}

}  // End namespace Foam

//****************************************************************************//


