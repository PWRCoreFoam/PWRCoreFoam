/*----------------------------------------------------------------------------*\
| Class:                                                                       |
|     Foam::rodMesh                                                            |
\*----------------------------------------------------------------------------*/

#ifndef rodMesh_H
#define rodMesh_H

#include "GeoMesh.H"
#include "subChannelMesh.H"


//****************************************************************************//

namespace Foam
{

/*----------------------------------------------------------------------------*\
|                       Class rodMesh Declaration                              |
\*----------------------------------------------------------------------------*/

class rodMesh
:
    public GeoMesh<subChannelMesh>
{

public:

    // Constructors

        // Construct from fvMesh
        explicit rodMesh(const subChannelMesh& mesh)
        :
            GeoMesh<subChannelMesh>(mesh)
        {}


    // Public member functions

        // Return size
        label size() const
        {
           return size(mesh_);
        }

        // Return size
        static label size(const Mesh& mesh)
        {
            return mesh.n_rods_sections();
        }

};

}  // End namespace Foam

//****************************************************************************//

#endif


