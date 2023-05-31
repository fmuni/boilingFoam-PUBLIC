scalar sumConduction = 0.0;
scalar sumConvection = 0.0;
scalar sumTot = 0.;

//- First, compute the heat fluxes in the fluid region
surfaceScalarField heatFluxD = -fvc::interpolate(kappa)*fvc::snGrad(T);

const surfaceScalarField::Boundary& patchHeatFluxD = heatFluxD.boundaryFieldRef();

if(mixture.found("writeQwall"))
{
    qwFluid[i] = -kappa*fvc::grad(T);
}


Info<< "\nWall heat fluxes [W] fluid region " << fluidRegions[i].name()
    << endl;
forAll(patchHeatFluxD, patchi)
{
    if(isA<processorFvPatchScalarField>(T.boundaryField()[patchi]))
    {
        continue;
    }

    scalar conduction = gSum
    (
        mesh.magSf().boundaryField()[patchi]
       *heatFluxD.boundaryField()[patchi]
    );

    scalar convection = 0.0;

    if(mesh.foundObject<surfaceScalarField>("phi"))
    {
        const surfaceScalarField& phi =
            mesh.lookupObject<surfaceScalarField>("phi");

        convection = gSum
        (
            rhoCp.boundaryField()[patchi]
            *T.boundaryField()[patchi]
            *phi.boundaryField()[patchi]
        );
    }


    Info<< mesh.boundary()[patchi].name()
        << " "
        << conduction
        << " "
        << convection
        << " "
        << conduction + convection
        << endl;

    sumConduction += conduction;
    sumConvection += convection;

}

sumTot = sumConduction + sumConvection;
Info<< "sum "
    << sumConduction
    << " "
    << sumConvection
    << " "
    << sumTot
    << nl << endl;

//- Compute maximum temperature
scalar maxTnew(gMax(T));

scalar Tres(mag(maxTFluid[i]-maxTnew)/(maxTFluid[i]+SMALL));

maxTFluid[i] = maxTnew;


Info<< "maxT = " << maxTFluid[i]  << "   residuals on T = "
    << Tres << nl << endl;

if(Pstream::master())
{
    std::ofstream outfile;

    outfile.open
    (
        fluidRegions[i].name() + "_heatFluxBalance.csv", std::ios_base::app
    );

    outfile << runTime.value() << ", "
            << sumConduction + sumConvection<< "\n";

    outfile.close();

    outfile.open
    (
        fluidRegions[i].name() + "_maxResT.csv", std::ios_base::app
    );

    outfile << runTime.value() << ", "
            << maxTFluid[i] << ", "
            << Tres << "\n";


}

if(mag(sumTot) < fluxTol && Tres < TTol)
{
    converged = true;
}
else
{
    converged = false;
}
