{

    kappa = alpha1*k1 + alpha2*k2;
    surfaceScalarField kappaf(fvc::interpolate(kappa));
    volScalarField energyExp(energySource*Tsat);

    #include "computeMassSource.H"



    fvScalarMatrix TEqn
     (
         fvm::ddt(rhoCp,T)
       + fvm::div(rhoCpPhi,T)
       - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T)
       - fvm::laplacian(kappaf,T)
       - fvm::Sp(energySource,T)
  //     - fvc::DDt(phi,p_rgh)
       ==
         - energyExp
     );


     TEqn.relax();
     TEqn.solve();

     Info<< "min/max(T) = " << min(T).value() << ", " << max(T).value() <<endl;

     //update continuity source with new T
     #include "computeMassSource.H"


}
