/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    chtCoalChemistryFoam

Description
    Solver for steady or transient fluid flow and solid heat conduction, with
    conjugate heat transfer between regions, compressible effects, turbulence,
    coal and limestone particle clouds, an energy source, and combustion.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermoMomentumTransportModel.H"
#include "psiReactionThermophysicalTransportModel.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "fixedGradientFvPatchFields.H" // added
#include "regionProperties.H" 		// added
#include "compressibleCourantNo.H" 	// added
#include "solidRegionDiffNo.H" 		// added
#include "solidThermo.H" 		// added

#include "basicThermoCloud.H"
#include "coalCloud.H"
#include "fvOptions.H"
#include "radiationModel.H"
#include "SLGThermo.H"
// #include "pimpleControl.H"
#include "coordinateSystem.H"
#include "pimpleMultiRegionControl.H"	// added
#include "pressureControl.H"
// #include "localEulerDdtScheme.H"
// #include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL						  // added
    #define CREATE_MESH createMeshesPostProcess.H		  // added
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    // #include "createMesh.H"
    #include "createMeshes.H"		// added
    #include "createFields.H"
    #include "initContinuityErrs.H"
    pimpleMultiRegionControl pimples(fluidRegions, solidRegions); // added
    #include "createFluidPressureControls.H" 			  // added
    // #include "createControl.H"
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"	// added
    // #include "createFieldRefs.H"
    #include "compressibleMultiRegionCourantNo.H"		   // added
    #include "solidRegionDiffusionNo.H"				   // added
    #include "setInitialMultiRegionDeltaT.H"			   // added

    //turbulence->validate();

    /*if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }*/

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    //while (pimple.run(runTime))
    while (pimples.run(runTime))
    {
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"		// added

 	#include "compressibleMultiRegionCourantNo.H"	// added
	#include "solidRegionDiffusionNo.H"		// added
	#include "setMultiRegionDeltaT.H"		// added

        /*if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }*/

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //rhoEffLagrangian = coalParcels.rhoEff() + limestoneParcels.rhoEff();
        //pDyn = 0.5*rho*magSqr(U);

        //coalParcels.evolve();

        //limestoneParcels.evolve();

        // #include "rhoEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        // while (pimple.loop())
        while (pimples.loop())
        {
           /*#include "UEqn.H"
            #include "YEqn.H"
            #include "EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
                thermophysicalTransport->correct();
            }*/

	    forAll(fluidRegions, i)
	    {
		Info<< "\nSolving for fluid region "
		    << fluidRegions[i].name() << endl;
		    #include "setRegionFluidFields.H"
		    #include "solveFluid.H"
	    }

	    forAll(solidRegions, i)
	    {
		Info<< "\nSolving for solid region "
		    << solidRegions[i].name() << endl;
		    #include "setRegionSolidFields.H"
		    #include "solveSolid.H"
	    }
        }

        // rho = thermo.rho();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
