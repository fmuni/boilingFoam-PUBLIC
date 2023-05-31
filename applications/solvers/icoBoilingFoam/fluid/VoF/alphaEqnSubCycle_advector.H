if (nOuterCorr > 1)
{
    // If nOuterCorrectors > 1 then for all but the first loop the advection
    // of alpha is done using an average, 0.5*phi+0.5*phiNew where phi is
    // the flux at the beginning of the time step and phiNew is the flux
    // estimate at the end of the time step from the previous outer
    // iteration. Similarly we use 0.5*U + 0.5*UNew in later iterations.
    if (oCorr == 0)
    {
        // To recalculate the alpha1 update in subsequent iterations, we
        // must store its current value before overwriting with the new
        // value
        alpha1.storePrevIter();
        // Storing initial phi and U for use in later outer iterations.
        phi.storePrevIter();
        U.storePrevIter();
    }
    else
    {
        // Resetting alpha1 to value before advection in first PIMPLE
        // iteration.
        alpha1 = alpha1.prevIter();

        // Temporarily setting U and phi with which to advect interface.
        U = 0.5*U.prevIter() + 0.5*U;
        phi = 0.5*phi.prevIter() + 0.5*phi;
    }
}

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
    
    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn_advector.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        rhoPhiCpSum += (runTime.deltaT()/totalDeltaT)*rhoCpPhi;
    }

    rhoPhi = rhoPhiSum;
    rhoCpPhi = rhoPhiCpSum;
}
else
{
    #include "alphaEqn_advector.H"
}

if (oCorr != 0)
{
    // Resetting U and phi to value at latest iteration.
    U = 2.0*U - U.prevIter();
    phi = 2.0*phi - phi.prevIter();
}

rho == alpha1*rho1 + alpha2*rho2;
rhoCp == alpha1*rho1*cp1 + alpha2*rho2*cp2;
