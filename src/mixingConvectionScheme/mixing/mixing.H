/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

Class
    Foam::mixing

Description
    mixing differencing scheme class.

SourceFiles
    mixing.C

\*---------------------------------------------------------------------------*/

#ifndef mixing_H
#define mixing_H

#include "limitedSurfaceInterpolationScheme.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mixing Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mixing
:
    public limitedSurfaceInterpolationScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const mixing&);


public:

    //- Runtime type information
    TypeName("mixing");


    // Constructors

        //- Construct from faceFlux
        mixing
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, faceFlux)
        {}

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        mixing
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, is)
        {}

        //- Construct from faceFlux and Istream
        mixing
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream&
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, faceFlux)
        {}


    // Member Functions

        //- Return the interpolation limiter
        virtual tmp<surfaceScalarField> limiter
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return tmp<surfaceScalarField>
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "mixingLimiter",
                        this->mesh().time().timeName(),
                        this->mesh()
                    ),
                    this->mesh(),
                    dimensionedScalar("mixingLimiter", dimless, 0.0)
                )
            );
        }

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights() const
        {

            //return pos(this->faceFlux_);

            return tmp<surfaceScalarField>
            (
              new surfaceScalarField
              (
                IOobject
                (
                  "weights",
                  this->mesh().time().timeName(),
                  this->mesh()
                ),
                this->mesh(),
                dimensionedScalar("weights", this->faceFlux_.dimensions(), 1.0)
              )
            );
        }

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return weights();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


