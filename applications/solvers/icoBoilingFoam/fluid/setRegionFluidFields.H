dynamicFvMesh& mesh = fluidRegions[i];

incompressibleTwoPhaseMixture& mixture(mixtureFluid[i]);

Boiling::interfaceProperties& interface(interfaceFluid[i]);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

volScalarField& rho = rhoFluid[i];
volVectorField& U = UFluid[i];
surfaceScalarField& phi = phiFluid[i];
surfaceScalarField& alphaPhi10 = alphaPhi10Fluid[i];

incompressible::turbulenceModel& turbulence = turbulenceFluid[i];
volScalarField& p_rgh = p_rghFluid[i];
volScalarField& p = pFluid[i];

volScalarField& T = TFluid[i];

volScalarField& kappa = kappaFluid[i];

volScalarField& mDot = mDotFluid[i];
volScalarField& mDot0 = mDot0Fluid[i];

const bool& alphaRestart = alphaRestartFluid[i];

const bool& useIsoAdvection = useIsoAdvectionFluid[i];

scalar& cumulativeContErr = cumulativeContErrFluid[i];

volScalarField& gh = ghFluid[i];
surfaceScalarField& ghf = ghfFluid[i];
const dimensionedScalar& ghRef = ghRefFluid[i];

fv::options& fvOptions = fluidFvOptions[i];

bool frozenFlow = frozenFlowFluid[i];

const label pRefCell = pRefCellFluid[i];
const scalar pRefValue = pRefValueFluid[i];

dictionary phase1 = mixture.subDict(mixture.phase1Name());
dictionary phase2 = mixture.subDict(mixture.phase2Name());
dictionary pcProperties = mixture.subDict("PhaseChangeProperties");

dimensionedScalar k1
(
    "kappa",
    dimensionSet(1,1,-3,-1,0,0,0),
    readScalar(phase1.lookup("kappa"))
);
dimensionedScalar cp1
(
    "Cp",
    dimensionSet(0,2,-2,-1,0,0,0),
    readScalar(phase1.lookup("Cp"))
);

dimensionedScalar k2
(
    "kappa",
    dimensionSet(1,1,-3,-1,0,0,0),
    readScalar(phase2.lookup("kappa"))
);
dimensionedScalar cp2
(
    "Cp",
    dimensionSet(0,2,-2,-1,0,0,0),
    readScalar(phase2.lookup("Cp"))
);

dimensionedScalar Tsat
(
    "Tsat",
    dimTemperature,
    readScalar(pcProperties.lookup("Tsat"))
);
dimensionedScalar R
(
    "R",
    dimensionSet(0,2,-2,-1,0,0,0),
    readScalar(pcProperties.lookup("R"))
);
dimensionedScalar sigmaEvap
(
    "sigmaEvap",
    dimless,
    readScalar(pcProperties.lookup("sigmaEvap"))
);

dimensionedScalar DmDot
(
    "DmDot",
    dimArea,
    readScalar(pcProperties.lookup("DmDot"))
);

dimensionedScalar hf
(
    "hf",
    dimensionSet(0,2,-2,0,0,0,0),
    readScalar(pcProperties.lookup("hf"))
);

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1*cp1 + alpha2*rho2*cp2,
    alpha1.boundaryField().types()
);
rhoCp.oldTime();


//Energy flux
surfaceScalarField rhoCpPhi
(
    IOobject
    (
        "rhoPhiCp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rhoCp)*phi
);

surfaceScalarField alphaPhiUn
(
    IOobject
    (
        "alphaPhiUn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(phi.dimensions(), Zero)
);

volScalarField continuitySource
(
    IOobject
    (
        "continuitySource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "continuitySource",
        dimensionSet(1,-3,-1,0,0,0,0),
        0.0
    )
);

volScalarField energySource
(
    IOobject
    (
        "energySource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "energySource",
        dimensionSet(1,-1,-3,-1,0,0,0),
        0.0
    )
);

dimensionedScalar evapRate("evapRate",dimensionSet(1,0,-1,0,0,0,0), 0.0);
dimensionedScalar condRate("condRate",dimensionSet(1,0,-1,0,0,0,0), 0.0);
dimensionedScalar totMdot("totMdot",dimensionSet(1,0,-1,0,0,0,0), 0.0);
dimensionedScalar mDotL("mDotL",dimensionSet(1,0,-1,0,0,0,0), 0.0);
dimensionedScalar mDotV("mDotV",dimensionSet(1,0,-1,0,0,0,0), 0.0);

pimpleControl& pimple(pimpleFluid[i]);

bool correctPhi
(
    pimple.dict().getOrDefault("correctPhi", mesh.dynamic())
);

bool checkMeshCourantNo
(
    pimple.dict().getOrDefault("checkMeshCourantNo", false)
);

bool moveMeshOuterCorrectors
(
    pimple.dict().getOrDefault("moveMeshOuterCorrectors", false)
);

IOMRFZoneList& MRF = MRFfluid[i];

surfaceScalarField& alphaPhi0(alphaPhi10Fluid[i]);

// MULES Correction
tmp<surfaceScalarField> talphaPhi1Corr0;

autoPtr<surfaceVectorField> Uf;

if (fluidRegions[i].dynamic())
{
    Info<< "Constructing face velocity Uf\n" << endl;

    Uf.reset
    (
        new surfaceVectorField
        (
            IOobject
            (
                "Uf",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(UFluid[i])
        )
    );
}

volScalarField rAU
(
    IOobject
    (
        "rAU",
        runTime.timeName(),
        fluidRegions[i],
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fluidRegions[i],
    dimensionedScalar("rAU", dimTime/dimDensity, 1)
);
