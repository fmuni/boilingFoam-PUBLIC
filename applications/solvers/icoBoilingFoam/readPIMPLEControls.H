    // We do not have a top-level mesh. Construct the fvSolution for
    // the runTime instead.
    fvSolution solutionDict(runTime);

    const dictionary& pimple = solutionDict.subDict("PIMPLE");

    const int nOuterCorr =
        pimple.getOrDefault<int>("nOuterCorrectors", 1);

    const int nMeshCorr =
        pimple.getOrDefault<int>("nMeshCorrectors", 1);

    bool converged(false);

    const scalar fluxTol(pimple.getOrDefault<double>("fluxTolerance", 1e-8));
    const scalar TTol(pimple.getOrDefault<double>("TTolerance", 1e-8));

    const Switch computeFluxes =
        pimple.getOrDefault<Switch>("computeFluxes", false);

    const Switch stopWhenConverged =
        pimple.getOrDefault<Switch>("stopWhenConverged", false);
