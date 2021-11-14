/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

\*---------------------------------------------------------------------------*/

#include "mixingConvectionScheme.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMatrices.H"

#include "subChannelMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
const surfaceInterpolationScheme<Type>&
mixingConvectionScheme<Type>::interpScheme() const
{
    return tinterpScheme_();
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
mixingConvectionScheme<Type>::interpolate
(
    const surfaceScalarField&,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return tinterpScheme_().interpolate(vf);
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
mixingConvectionScheme<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return faceFlux*interpolate(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type> >
mixingConvectionScheme<Type>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<surfaceScalarField> tweights = tinterpScheme_().weights(vf);
    const surfaceScalarField& weights = tweights();

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            faceFlux.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm();


    const Foam::subChannelMesh & mesh_ =
      faceFlux.db().lookupObject<Foam::subChannelMesh>("data");

    labelList gap_list_ = mesh_.get_gapList();

    const labelUList& l = fvm.lduAddr().lowerAddr();
    const labelUList& u = fvm.lduAddr().upperAddr();

    forAll(fvm.lower(), faceI)
    {
      fvm.lower()[faceI] = 0.0;
      fvm.upper()[faceI] = 0.0;

      fvm.diag()[l[faceI]] = 0.0;
      fvm.diag()[u[faceI]] = 0.0;
    }

    forAll(gap_list_, gap_i)
    {
      label faceI = gap_list_[gap_i];

      fvm.lower()[faceI] = -1.0 * faceFlux[faceI];
      fvm.upper()[faceI] = -1.0 * faceFlux[faceI];

      fvm.diag()[l[faceI]] += faceFlux[faceI];
      fvm.diag()[u[faceI]] += faceFlux[faceI];
    }


    forAll(vf.boundaryField(), patchI)
    {
        const fvPatchField<Type>& psf = vf.boundaryField()[patchI];
        const fvsPatchScalarField& patchFlux = faceFlux.boundaryField()[patchI];
        const fvsPatchScalarField& pw = weights.boundaryField()[patchI];

        //fvm.internalCoeffs()[patchI] = patchFlux;
        //fvm.boundaryCoeffs()[patchI] = patchFlux;
        fvm.internalCoeffs()[patchI] = patchFlux*psf.valueInternalCoeffs(pw);
        fvm.boundaryCoeffs()[patchI] = -patchFlux*psf.valueBoundaryCoeffs(pw);
    }

    if (tinterpScheme_().corrected())
    {
        fvm += fvc::surfaceIntegrate(faceFlux*tinterpScheme_().correction(vf));
    }

    return tfvm;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
mixingConvectionScheme<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tConvection
    (
        fvc::surfaceIntegrate(flux(faceFlux, vf))
    );

    tConvection().rename
    (
        "convection(" + faceFlux.name() + ',' + vf.name() + ')'
    );

    return tConvection;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //


