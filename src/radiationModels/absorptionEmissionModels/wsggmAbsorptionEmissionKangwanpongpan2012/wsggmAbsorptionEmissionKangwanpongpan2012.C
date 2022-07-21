/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "wsggmAbsorptionEmissionKangwanpongpan2012.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "basicSpecieMixture.H"
#include "psiReactionThermo.H"
#include "surfaceFields.H"
#include "symmetryFvPatch.H"
#include "cyclicFvPatch.H"
#include "processorFvPatch.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiationModels
    {
	namespace absorptionEmissionModels
	{
        	defineTypeNameAndDebug(wsggmAbsorptionEmissionKangwanpongpan2012, 0);

	        addToRunTimeSelectionTable
	        (
	            absorptionEmissionModel,
	            wsggmAbsorptionEmissionKangwanpongpan2012,
	            dictionary
	        );
	}
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::wsggmAbsorptionEmissionKangwanpongpan2012::wsggmAbsorptionEmissionKangwanpongpan2012
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& modelName
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    thermo_(mesh.lookupObject<fluidThermo>(basicThermo::dictName)),
    pathLength_(coeffsDict_.lookup("pathLength")),
    meanBeamPathAutoCalcMode_(coeffsDict_.lookupOrDefault<bool>("meanBeamPathAutoCalcMode", false)),
    sector_(coeffsDict_.lookup("sector"))
{

    if (!isA<basicSpecieMixture>(thermo_))
    {
        FatalErrorInFunction
            << "Model requires a multi-component thermo package"
            << abort(FatalError);
    }

    label nD = mesh.nGeometricD();

    if (nD == 3)
    {

        if (meanBeamPathAutoCalcMode_)
        {
            scalar totVolume = gSum(mesh.V());

            scalar totArea = 0.0;
            forAll(mesh.boundary(),patchI)
            {
                if ( (!isA<processorFvPatch>(mesh.boundary()[patchI]))
                    && (!isA<symmetryFvPatch>(mesh.boundary()[patchI]))
                    && (!isA<cyclicFvPatch>(mesh.boundary()[patchI])) )
                  {
	                    totArea += sum(mesh.magSf().boundaryField()[patchI]);
	              }
            }
            reduce(totArea, sumOp<scalar>());

            pathLength_.value() = 3.6*(360.0/sector_.value())*totVolume/totArea;
            Info << "using the computed pathLength: " <<  pathLength_ << endl;
        }
        else
        {
            Info << "using the user-provided pathLength: "
                 <<  pathLength_ << endl;
        }

    }
    else if (nD == 2)
    {
        WarningInFunction
            << "Case is 2D, pathLength is not strictly applicable, \n using the User provided  pathLength: "
            <<  pathLength_
            << endl;
    }
    else
    {
        FatalErrorInFunction
            << "Case is not 3D or 2D, pathLength is not applicable"
            << exit(FatalError);
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::wsggmAbsorptionEmissionKangwanpongpan2012::~wsggmAbsorptionEmissionKangwanpongpan2012()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::wsggmAbsorptionEmissionKangwanpongpan2012::aCont(const label bandI) const
{
    const basicSpecieMixture& mixture =
        dynamic_cast<const basicSpecieMixture&>(thermo_);

    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();

    label indexCO2 = mixture.species()["CO2"];
    label indexH2O = mixture.species()["H2O"];

    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "aCont" + name(bandI),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("a", dimless/dimLength, 0.0)
        )
    );

    volScalarField& a = ta.ref();

    scalar Tref = 2000;

    //- See Table 3 from reference
    double CK_[4][3] = {
    {	0.0429,	  0.0093, -0.0018 },
    {	0.3647,	  0.0790, -0.0150 },
    {	3.7144,   0.2565, -0.0509 },
    { 105.3100,	-39.2650,  6.0877 }
    };

    //- See Table 4 from reference
    double C_[24][3] = {
    {  0.3947,  -0.1214,  0.0243 },
    { -0.4512,   1.1420, -0.2296 },
    {  0.1492,  -5.2222,  1.0115 },
    {  1.8824,   9.1820, -1.7493 },
    { -2.3284,  -6.9298,  1.3038 },
    {  0.7698,   1.9063, -0.3549 },
    { -0.4974,   0.1092, -0.0179 },
    {  6.8986,  -2.3198,  0.4077 },
    {-19.9880,   8.0021, -1.4482 },
    { 26.2080, -11.0070,  2.0311 },
    {-16.4400,   7.1199, -1.3278 },
    {  3.9847,  -1.7876,  0.3349 },
    {  0.3189,  -0.0720,  0.0158 },
    { -0.7222,   1.0304, -0.2478 },
    {  1.5053,  -1.9350,  0.5931 },
    { -1.8378,   1.6332, -0.6619 },
    {  1.0337,  -0.7798,  0.3857 },
    { -0.2107,   0.1782, -0.0933 },
    {  0.1648,   0.0329, -0.0095 },
    { -0.6012,   0.6942, -0.0687 },
    {  2.0308,  -3.0960,  0.3691 },
    { -3.4361,   4.7494, -0.5919 },
    {  2.5803,  -3.1714,  0.4017 },
    { -0.7069,   0.7869, -0.1003 }
    };

    forAll(a, celli)
    {
	    a[celli]=0.0;
	    scalar invWt = 0.0;
            forAll(mixture.Y(), s)
            {
                invWt += mixture.Y(s)[celli]/mixture.Wi(s);
            }

	    scalar XkCO2 = mixture.Y(indexCO2)[celli]/(mixture.Wi(indexCO2)*invWt);
	    scalar XkH2O = mixture.Y(indexH2O)[celli]/(mixture.Wi(indexH2O)*invWt);

	    scalar mr = XkH2O/XkCO2;

	    double emissivityCoeffs_[4];
	    double fittingFactors_[4][6];

	    //- Construct 1D array of absorption coeffs. (as function of molar ratio)
	    for(int i=0; i<4; i++)
	    {
	        emissivityCoeffs_[i] = CK_[i][0] + CK_[i][1]*pow(mr, 1) + CK_[i][2]*pow(mr, 2);
	    }

	    //- Construct 2D array of polynomial coeffs. (as function of molar ratio)
	    for(int i=0; i<4; i++)
	    {
		for(int j=0; j<6; j++)
		{
		    fittingFactors_[i][j] = C_[6*i+j][0] + C_[6*i+j][1]*pow(mr, 1) + C_[6*i+j][2]*pow(mr, 2);
		}
	    }

 	    scalar pressurePathLength = (XkCO2*paToAtm(p[celli]) + XkH2O*paToAtm(p[celli]))*pathLength_.value();

	    //- Data shown in charts is valid up to T = 2500K
	    scalar limitedDimlessTemperature = min(T[celli]/Tref, 1.25);

	    double emissivity = 0.0;
	    for(int i=0; i<4; i++)
            {
  	        double weightingFactor = 0.0;

		for(int j=0; j<6; j++)
		{
		    weightingFactor += fittingFactors_[i][j]*pow(limitedDimlessTemperature, (j));
		}
		emissivity += weightingFactor * (1 - exp( (-1) *emissivityCoeffs_[i] * pressurePathLength));
            }

	    emissivity = min(emissivity,0.9999);
	    a[celli] = (-1)* log(1-emissivity) / pathLength_.value();
    }

    ta.ref().correctBoundaryConditions(); //correct BC's
    return ta;
}

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::wsggmAbsorptionEmissionKangwanpongpan2012::eCont(const label bandI) const
{
    return aCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::wsggmAbsorptionEmissionKangwanpongpan2012::ECont(const label bandI) const
{
    tmp<volScalarField> tE
    (
        new volScalarField
        (
            IOobject
            (
                "ECont" + name(bandI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );

    return tE;

}


// ************************************************************************* //
