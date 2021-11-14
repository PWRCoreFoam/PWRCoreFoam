/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

Class
    Foam::fv::mixingConvectionScheme

Description
    Basic second-order convection using face-gradients and Gauss' theorem.

SourceFiles
    mixingConvectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef mixingConvectionScheme_H
#define mixingConvectionScheme_H

#include "convectionScheme.H"

#include "subChannelMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

//- Temporary debug switch to provide warning about backward-compatibility
//  issue with setting div schemes for steady-state
extern int warnUnboundedGaussMixing;

/*---------------------------------------------------------------------------*\
                       Class mixingConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mixingConvectionScheme
:
    public fv::convectionScheme<Type>
{
    // Private data

        tmp<surfaceInterpolationScheme<Type> > tinterpScheme_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        mixingConvectionScheme(const mixingConvectionScheme&);

        //- Disallow default bitwise assignment
        void operator=(const mixingConvectionScheme&);


public:

    //- Runtime type information
    TypeName("Mixing");


    // Constructors

        //- Construct from flux and interpolation scheme
        mixingConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const tmp<surfaceInterpolationScheme<Type> >& scheme
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            tinterpScheme_(scheme)
        {}

        //- Construct from flux and Istream
        mixingConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            tinterpScheme_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            )
        {
            is.rewind();
            word bounded(is);

            if
            (
                warnUnboundedGaussMixing
             && word(mesh.ddtScheme("default")) == "steadyState"
             && bounded != "bounded"
            )
            {
                fileNameList controlDictFiles(findEtcFiles("controlDict"));

                IOWarningIn("mixingConvectionScheme", is)
                    << "Unbounded 'Mixing' div scheme used in "
                       "steady-state solver, use 'bounded Gauss' "
                       "to ensure boundedness.\n"
                    << "    To remove this warning switch off "
                    << "'boundedGauss' in "
                    << controlDictFiles[controlDictFiles.size()-1]
                    << endl;
            }
        }


    // Member Functions

        const surfaceInterpolationScheme<Type>& interpScheme() const;

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > interpolate
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > flux
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<fvMatrix<Type> > fvmDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mixingConvectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


