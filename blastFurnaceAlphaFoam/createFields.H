#include "createRDeltaT.H"

#include "readGravitationalAcceleration.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

SLGThermo slgThermo(mesh, thermo);

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species()
        << exit(FatalIOError);
}

volScalarField& p = thermo.p();
const volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField rhoc
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

// lagrangian effective density field - used externally (optional)
volScalarField rhoEffLagrangian
(
    IOobject
    (
        "rhoEffLagrangian",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimDensity, 0)
);

// dynamic pressure field - used externally (optional)
volScalarField pDyn
(
    IOobject
    (
        "pDyn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimPressure, 0)
);

// radiative source term
volScalarField radiativeSourceTerm
(
    IOobject
    (
        "radiativeSourceTerm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
   mesh,
   dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
);

// (1) rename field U to Uc
Info<< "\nReading field U\n" << endl;
volVectorField Uc
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// (2) Comment: #include "compressibleCreatePhi.H"
// (3) add flux field phic
Info<< "Reading/calculating face flux field phic\n" << endl;

surfaceScalarField phic
(
    IOobject
    (
	"phi",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    linearInterpolate(Uc) & mesh.Sf()
);

mesh.setFluxRequired(p.name());

// (4) add volume fraction field alphac
volScalarField alphac
(
    IOobject
    (
	"alpha",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, Zero)
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

// (5) Replace K=K(U) for Kc=Kc(Uc)
Info<< "Creating field kinetic energy Kc\n" << endl;
volScalarField Kc("K", 0.5*magSqr(Uc));

if (Uc.nOldTimes())
{
    volVectorField* Ucold = &Uc.oldTime();
    volScalarField* Kcold = &Kc.oldTime();
    *Kcold == 0.5*magSqr(*Ucold);

    while (Ucold->nOldTimes())
    {
	Ucold = &Ucold->oldTime();
	Kcold = &Kcold->oldTime();
        *Kcold == 0.5*magSqr(*Ucold);
    }
}

#include "createMRF.H"
#include "createClouds.H"

// (6) calculate minimum volume fraction alphacMin
scalar alphacMin
(
    1.0
  - readScalar
    (
	biomassMPPICParcels.particleProperties().subDict("constantProperties")
	.lookup("alphaMax")
    )
  - readScalar
    (
	ironOreMPPICParcels.particleProperties().subDict("constantProperties")
	.lookup("alphaMax")
    )
  - readScalar
    (
	limestoneMPPICParcels.particleProperties().subDict("constantProperties")
	.lookup("alphaMax")
    )
);

// (7) update alphac from the particle locations
alphac = max(1.0 - biomassMPPICParcels.theta() - ironOreMPPICParcels.theta()
	   - limestoneMPPICParcels.theta(), alphacMin);
alphac.correctBoundaryConditions();

surfaceScalarField alphacf("alphaf", fvc::interpolate(alphac));
Info<<"alphacf: "<<alphacf<<endl;
surfaceScalarField rhocPhic
(
    IOobject
    (
	"rhoPhi",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    fvc::interpolate(rhoc)*phic
);

surfaceScalarField alphaRhoPhic
(
    IOobject
    (
	"alphaRhoPhi",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    alphacf*rhocPhic
);

Info<< "Creating turbulence model\n" << endl;
autoPtr<blastFurnaceAlphaFoamMomentumTransportModel> turbulence
(
    blastFurnaceAlphaFoamMomentumTransportModel::New
    (
	alphac,
        rhoc,
        Uc,
	alphaRhoPhic,
        phic,
        thermo
    )
);

Info<< "Creating thermophysical transport model\n" << endl;
//autoPtr<psiReactionThermophysicalTransportModel> thermophysicalTransport
autoPtr<blastFurnaceAlphaFoamPsiReactionThermophysicalTransportModel> thermophysicalTransport
(
    //psiReactionThermophysicalTransportModel::New(turbulence(), thermo)
    blastFurnaceAlphaFoamPsiReactionThermophysicalTransportModel::New(turbulence(), thermo)
);

Info<< "Creating combustion model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> combustion
(
    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
);

#include "createRadiationModel.H"
#include "createFvOptions.H"
