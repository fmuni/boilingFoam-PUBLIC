if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );

    surfaceScalarField rhoPhiCpSum("rhoPhiCpSum", 0.0*rhoCpPhi);

    tmp<volScalarField> trSubDeltaT;

    // if (LTS)
    // {
    //     trSubDeltaT =
    //         fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
    // }

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        rhoPhiCpSum += (runTime.deltaT()/totalDeltaT)*rhoCpPhi;
    }

    rhoPhi = rhoPhiSum;
    rhoCpPhi = rhoPhiCpSum;
}
else
{
    #include "alphaEqn.H"
}

rho == alpha1*rho1 + alpha2*rho2;
rhoCp == alpha1*rho1*cp1 + alpha2*rho2*cp2;
