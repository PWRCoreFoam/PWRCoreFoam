/*----------------------------------------------------------------------------*\
| Class:                                                                       |
|     Foam::GeometricRodField                                                  |
\*----------------------------------------------------------------------------*/

#ifndef GeometricRodField_H
#define GeometricRodField_H

#include "DimensionedField.H"
#include "dimensionedTypes.H"
#include "FieldField.H"
#include "lduInterfaceFieldPtrsList.H"
#include "LduInterfaceFieldPtrsList.H"
#include "pointFieldFwd.H"
#include "regIOobject.H"


//****************************************************************************//

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type, class GeoMesh>
class GeometricRodField;


/*----------------------------------------------------------------------------*\
|                       Class GeometricRodField Declaration                    |
\*----------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
class GeometricRodField
:
    public DimensionedField<Type, GeoMesh>
{

public:

    // Public typedefs

        typedef typename GeoMesh::Mesh Mesh;
        typedef typename GeoMesh::BoundaryMesh BoundaryMesh;

        typedef DimensionedField<Type, GeoMesh> DimensionedInternalField;
        typedef Field<Type> InternalField;


private:

    // Private data

        // Current time index used to trigger the storing of the old-time values
        mutable label timeIndex_;

        // Pointer to old time field
        mutable GeometricRodField<Type, GeoMesh>* field0Ptr_;

        // Pointer to previous iteration (used for under-relaxation)
        mutable GeometricRodField<Type, GeoMesh>* fieldPrevIterPtr_;

        const Mesh& subChannelMesh_;

        const labelList & rodEdgesList_;
        const labelList & rodEdgesListOnAssemblyEdge_;
        const labelList & rodEdgesListOnCoreEdge_;

        const scalarList & rodEdgesLengthList_;

        const pointField & rodEdgesCenterList_;

        scalar z_min;
        scalar z_max;


    // Private member functions

        void calc_z_min_max();


public:

    // Runtime type information

        TypeName("GeometricRodField");


    // Public typedefs

        typedef typename Field<Type>::cmptType cmptType;


    // Static member functions

        // Return a null geometric field
        inline static const GeometricRodField<Type, GeoMesh>& null();


    // Constructors

        // Constructor given IOobject, mesh, dimensioned<Type>.
        GeometricRodField
        (
            const IOobject&,
            const Mesh&,
            const dimensioned<Type>&
        );

        // Constructor given IOobject and mesh.
        GeometricRodField
        (
            const IOobject&,
            const Mesh&
        );


    // Destructor

    virtual ~GeometricRodField();


    // Public member functions

        void writeRodField() const;
        void writeRodOnAssemblyEdge() const;
        void writeRodOnCoreEdge() const;

        void writeRodField_Test() const;

        inline const Mesh& rodMesh() const;

        inline const labelList & get_rodEdgesList() const;

        inline const scalarList & get_rodEdgesLengthList() const;

        inline const pointField & get_rodEdgesCenterList() const;

        inline scalar get_z_min() const;

        inline scalar get_z_max() const;


    // Member operators


    // IOstream operators


};

}  // End namespace Foam

//****************************************************************************//

#include "GeometricRodFieldI.H"

#ifdef NoRepository
#   include "GeometricRodField.C"
#endif

//****************************************************************************//

#endif


