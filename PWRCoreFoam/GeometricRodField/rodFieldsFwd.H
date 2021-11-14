/*----------------------------------------------------------------------------*\
| Class:                                                                       |
|     Foam::rodField                                                           |
\*----------------------------------------------------------------------------*/

#ifndef rodFieldsFwd_H
#define rodFieldsFwd_H

#include "fieldTypes.H"
#include "GeometricRodField.H"
#include "rodMesh.H"


//****************************************************************************//

namespace Foam
{

class rodMesh;

template<class Type, class GeoMesh>
class GeometricRodField;

typedef GeometricRodField<scalar, rodMesh> rodScalarField;

}  // End namespace Foam

//****************************************************************************//

#endif


