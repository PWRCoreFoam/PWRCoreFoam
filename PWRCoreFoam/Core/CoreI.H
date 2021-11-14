/*----------------------------------------------------------------------------*\
| Class:                                                                       |
|     Foam::Core                                                               |
\*----------------------------------------------------------------------------*/


//****************************************************************************//

namespace Foam
{

/*----------------------------------------------------------------------------*\
|                       Member functions                                       |
\*----------------------------------------------------------------------------*/

inline const Foam::IOdictionary & Foam::Core::get_core_dict()
{
  return coreDict_;
}


inline const Foam::IOdictionary & Foam::Core::get_core_dict() const
{
  return coreDict_;
}


inline const Foam::subChannelMesh & Foam::Core::mesh()
{
  return mesh_;
}


inline const Foam::subChannelMesh & Foam::Core::mesh() const
{
  return mesh_;
}


inline Foam::volVectorField & Foam::Core::get_momentum_source()
{
  return momentumSource_;
}


inline const Foam::volVectorField & Foam::Core::get_momentum_source() const
{
  return momentumSource_;
}


inline Foam::volScalarField & Foam::Core::get_energy_source()
{
  return energySource_;
}


inline const Foam::volScalarField & Foam::Core::get_energy_source() const
{
  return energySource_;
}


inline Foam::volScalarField & Foam::Core::get_Re()
{
  return Re_;
}


inline const Foam::volScalarField & Foam::Core::get_Re() const
{
  return Re_;
}


inline Foam::surfaceScalarField & Foam::Core::get_phi_mix()
{
  return phi_mix_;
}


inline const Foam::surfaceScalarField & Foam::Core::get_phi_mix() const
{
  return phi_mix_;
}

}  // End namespace Foam

//****************************************************************************//


