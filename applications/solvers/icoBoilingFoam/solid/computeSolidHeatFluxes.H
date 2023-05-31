scalar sumConduction = 0.0;

//- First, compute the heat fluxes in the fluid region
surfaceScalarField heatFluxD =
    -fvc::interpolate(mag(thermo.Kappa()))*fvc::snGrad(thermo.T());

const surfaceScalarField::Boundary& patchHeatFluxD = heatFluxD.boundaryFieldRef();


Info<< "\nWall heat fluxes [W] solid region " << solidRegions[i].name()
    << endl;
forAll(patchHeatFluxD, patchi)
{
    if(isA<processorFvPatchScalarField>(thermo.T().boundaryField()[patchi]))
    {
        continue;
    }

    scalar conduction = gSum
    (
        mesh.magSf().boundaryField()[patchi]
       *heatFluxD.boundaryField()[patchi]
    );



    Info<< mesh.boundary()[patchi].name()
        << " "
        << conduction
        << endl;

    sumConduction += conduction;

}

Info<< "sum "
    << sumConduction
    << nl << endl;

//- Compute maximum temperature
scalar maxTnew(gMax(thermo.T()));

scalar Tres(mag(maxTSolid[i]-maxTnew)/(maxTSolid[i]+SMALL));

maxTSolid[i] = maxTnew;

Info<< "maxT = " << maxTSolid[i] << "   residuals on T = " << Tres << endl;

if(Pstream::master())
{
    std::ofstream outfile;

    outfile.open
    (
        solidRegions[i].name() + "_heatFluxBalance.csv", std::ios_base::app
    );

    outfile << runTime.value() << ", "
            << sumConduction << "\n";

    outfile.close();

    outfile.open
    (
        solidRegions[i].name() + "_maxResT.csv", std::ios_base::app
    );

    outfile << runTime.value() << ", "
            << maxTSolid[i] << ", "
            << Tres << "\n";


}

if(mag(sumConduction) < fluxTol && Tres < TTol)
{
    converged = true;
}
else
{
    converged = false;
}
