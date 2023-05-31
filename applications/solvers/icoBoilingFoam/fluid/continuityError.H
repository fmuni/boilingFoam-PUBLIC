
volScalarField contErr(fvc::ddt(rho) + fvc::div(rhoPhi));

//volScalarField contErr(fvc::div(phi) - vDotCond + vDotEvap);

scalar sumLocalContErr = runTime.deltaTValue()*
    mag(contErr)().weightedAverage(mesh.V()).value();

scalar globalContErr = runTime.deltaTValue()*
    contErr.weightedAverage(mesh.V()).value();

cumulativeContErr += globalContErr;

Info<< "time step continuity errors : sum local = " << sumLocalContErr
    << ", global = " << globalContErr
    << ", cumulative = " << cumulativeContErr
    << endl;
