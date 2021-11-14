/*----------------------------------------------------------------------------*\
| Application                                                                  |
\*----------------------------------------------------------------------------*/

#include "Core.H"
#include "fvCFD.H"
#include "rhoThermo.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "IOporosityModelList.H"
#include "fixedFluxPressureFvPatchScalarField.H"

#include "rodFieldsFwd.H"
#include "subChannelMesh.H"


//****************************************************************************//

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    #include "createSubChannelMesh.H"

    #include "readGravitationalAcceleration.H"
    #include "createFields.H"

    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"

    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    pimpleControl pimple(mesh);

    //************************************************************************//

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"

            #include "EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        #include "writeData.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


//****************************************************************************//


