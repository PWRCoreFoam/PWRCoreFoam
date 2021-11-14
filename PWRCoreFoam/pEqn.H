{

    rho = thermo.rho();

    volScalarField rAU(1.0/UEqn().A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn().H();


    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(rho*HbyA) & mesh.Sf())
      + (fvc::interpolate(rho*momentumSource*rAU) & mesh.Sf())
    );

// Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        // Pressure corrector
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rhorAUf, p) == fvc::div(phiHbyA)
        );
    
        pEqn.setReference(pRefCell, pRefValue);
    
        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
    
        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

#include "compressibleContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
    p.relax();

    // Here maybe we can use fvc::reconstruct
    U = HbyA - rAU*fvc::grad(p) + rAU*momentumSource;

    U.correctBoundaryConditions();

    fvOptions.correct(U);

    K = 0.5*magSqr(U);

    if (thermo.dpdt())
    {
      dpdt = fvc::ddt(p);
    }

}


