/*----------------------------------------------------------------------------*\
| Class:                                                                       |
|     Foam::Core                                                               |
\*----------------------------------------------------------------------------*/

#ifndef Core_H
#define Core_H

#include "Fuel.H"
#include "fvCFD.H"
#include "IOdictionary.H"
#include "rhoThermo.H"
#include "subChannelMesh.H"


//****************************************************************************//

namespace Foam
{

/*----------------------------------------------------------------------------*\
|                       Class Core Declaration                                 |
\*----------------------------------------------------------------------------*/

class Core
:
    public Fuel
{

public:

    // Public typedefs


private:

    // Private data

        const IOdictionary coreDict_;

        const subChannelMesh & mesh_;

        const volVectorField & U_;

        const rhoThermo & thermo_;

        volVectorField momentumSource_;

        volScalarField energySource_;

        volScalarField Re_;

        surfaceScalarField phi_mix_;

        scalar beta_;

        volScalarField * test_variable_1;


    // Private member functions

        void calc_Re();

        void update_momentum_source();

        void calc_phi_mix();

        void update_energy_source();

        volScalarField * test_func_1() const;


public:

    // Runtime type information

        TypeName("Core");


    // Public typedefs


    // Static member functions


    // Constructors

        // Construct from subChannelMesh
        Core
        (
          const subChannelMesh & mesh,
          const volVectorField & U,
          const rhoThermo & thermo
        );

    // Destructor


    // Public member functions

        void update();

        void write();

        inline const IOdictionary & get_core_dict();
        inline const IOdictionary & get_core_dict() const;

        inline const subChannelMesh & mesh();
        inline const subChannelMesh & mesh() const;

        inline volVectorField & get_momentum_source();
        inline const volVectorField & get_momentum_source() const;

        inline volScalarField & get_energy_source();
        inline const volScalarField & get_energy_source() const;

        inline volScalarField & get_Re();
        inline const volScalarField & get_Re() const;

        inline surfaceScalarField & get_phi_mix();
        inline const surfaceScalarField & get_phi_mix() const;


    // Member operators


    // IOstream operators


};

}  // End namespace Foam

//****************************************************************************//

#include "CoreI.H"

//****************************************************************************//

#endif


