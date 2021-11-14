/*----------------------------------------------------------------------------*\
| Class:                                                                       |
|     Foam::subChannelMesh                                                     |
\*----------------------------------------------------------------------------*/


//****************************************************************************//

namespace Foam
{

/*----------------------------------------------------------------------------*\
|                       Member functions                                       |
\*----------------------------------------------------------------------------*/

inline label
subChannelMesh::coreAssemblyDimension()
{
  return CoreMapDimension_;
}


inline label
subChannelMesh::coreAssemblyDimension() const
{
  return CoreMapDimension_;
}


inline label
subChannelMesh::assemblyRodsDimension()
{
  return AssemblyRodsDimension_;
}


inline label
subChannelMesh::assemblyRodsDimension() const
{
  return AssemblyRodsDimension_;
}


inline scalar
subChannelMesh::assemblyWidth()
{
  return AssemblyWidth_;
}


inline scalar
subChannelMesh::assemblyWidth() const
{
  return AssemblyWidth_;
}


inline label
subChannelMesh::centerAssemblyXIndex()
{
  return CenterAssemblyXIndex_;
}


inline label
subChannelMesh::centerAssemblyXIndex() const
{
  return CenterAssemblyXIndex_;
}


inline label
subChannelMesh::centerAssemblyYIndex()
{
  return CenterAssemblyYIndex_;
}


inline label
subChannelMesh::centerAssemblyYIndex() const
{
  return CenterAssemblyYIndex_;
}


inline label
subChannelMesh::n_rods_sections()
{
  if (!rodEdgesListPtr_)
  {
    FatalErrorIn("void Foam::subChannelMesh::n_rods_sections()")
              << "Rods have not been built !"
              << exit(FatalError);
  }

  return (*rodEdgesListPtr_).size();
}


inline label
subChannelMesh::n_rods_sections() const
{
  if (!rodEdgesListPtr_)
  {
    FatalErrorIn("void Foam::subChannelMesh::n_rods_sections()")
              << "Rods have not been built !"
              << exit(FatalError);
  }

  return (*rodEdgesListPtr_).size();
}


inline scalarList &
subChannelMesh::get_flow_area() const
{
  return flow_area_;
}


inline scalarList &
subChannelMesh::get_hydraulic_diameter() const
{
  return hydraulic_diameter_;
}


inline const labelListList &
subChannelMesh::get_rod_assembly_index_x_y() const
{
  return (*rod_assembly_index_X_Y_);
}


inline const labelListList &
subChannelMesh::get_rod_rod_index_x_y() const
{
  return (*rod_rod_index_X_Y_);
}


inline const scalarListList &
subChannelMesh::get_core_assembly_power_distribution_factor() const
{
  return core_assembly_power_distribution_;
}


inline const scalarListList &
subChannelMesh::get_rod_pin_power_distribution_factor() const
{
  return rod_pin_power_distribution_factor_;
}


inline labelList &
subChannelMesh::get_core_type() const
{
  return (*core_type_);
}


inline scalar
subChannelMesh::mixing_coefficient()
{
  return mixing_coefficient_;
}


inline scalar
subChannelMesh::mixing_coefficient() const
{
  return mixing_coefficient_;
}


inline scalar
subChannelMesh::get_pitch()
{
  return pitch_;
}


inline scalar
subChannelMesh::get_pitch() const
{
  return pitch_;
}

}  // End namespace Foam

//****************************************************************************//


