/*----------------------------------------------------------------------------*\
| Class:                                                                       |
|     Foam::GeometricRodField                                                  |
\*----------------------------------------------------------------------------*/


//****************************************************************************//

namespace Foam
{

/*----------------------------------------------------------------------------*\
|                       Member functions                                       |
\*----------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
inline const typename GeoMesh::Mesh&
GeometricRodField<Type, GeoMesh>::rodMesh() const
{
    return subChannelMesh_;
}


template<class Type, class GeoMesh>
inline const labelList &
GeometricRodField<Type, GeoMesh>::get_rodEdgesList() const
{
  return rodEdgesList_;
}


template<class Type, class GeoMesh>
inline const scalarList &
GeometricRodField<Type, GeoMesh>::get_rodEdgesLengthList() const
{
  return rodEdgesLengthList_;
}


template<class Type, class GeoMesh>
inline const pointField &
GeometricRodField<Type, GeoMesh>::get_rodEdgesCenterList() const
{
  return rodEdgesCenterList_;
}


template<class Type, class GeoMesh>
inline scalar
GeometricRodField<Type, GeoMesh>::get_z_min() const
{
  return z_min;
}


template<class Type, class GeoMesh>
inline scalar
GeometricRodField<Type, GeoMesh>::get_z_max() const
{
  return z_max;
}

}  // End namespace Foam

//****************************************************************************//


