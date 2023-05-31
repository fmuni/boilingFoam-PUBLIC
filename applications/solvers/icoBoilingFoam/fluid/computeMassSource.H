{
   //alpha cutoff for source term distribution
   //no source term in cells where cutoff < alpha1 < 1 - cutoff

   scalar cutoff = 1.0e-3;

   dimensionedScalar T0("0", dimTemperature, 0);

   volScalarField intArea = mag(fvc::grad(alpha1));

   dimensionedScalar intAreaTot = fvc::domainIntegrate(intArea);
   dimensionedScalar intAreaTotLiquid = fvc::domainIntegrate(intArea*alpha1);
   dimensionedScalar intAreaTotVap = fvc::domainIntegrate(intArea*alpha2);

   dimensionedScalar Nl("Nl", dimensionSet(0,0,0,0,0,0,0), 2.0);
   dimensionedScalar Nv("Nv", dimensionSet(0,0,0,0,0,0,0), 2.0);

   if(intAreaTotLiquid.value() != 0) { Nl = intAreaTot / intAreaTotLiquid; }
   if(intAreaTotVap.value() != 0) { Nv = intAreaTot / intAreaTotVap; }

   dimensionedScalar RintEvap =
    ((2.0*sigmaEvap)/(2.0-sigmaEvap)) * (rho2*hf)
    / pow(2.0*M_PI*R*pow(Tsat,3),0.5);

   //Non smoothed evaporation rate [Kg/m^3 s]
   mDot0  = Nl*alpha1*RintEvap*intArea*(T - Tsat);
   energySource =  - Nl*RintEvap*hf*intArea*alpha1;
   energyExp = energySource*Tsat;

   #include "microlayerMassSource.H"

   //Total evaporation rate [Kg/s]
   evapRate = fvc::domainIntegrate(mDot0);

   //Smooth mDot
   fvScalarMatrix mDotEqn
   (
      fvm::Sp(scalar(1),mDot) - fvm::laplacian(DmDot,mDot) == mDot0
   );

   mDotEqn.solve();

   //Total evaporation rate after smoothing, should be equal to evapRate
   totMdot = fvc::domainIntegrate(mDot);

   dimensionedScalar intMdotVapour
   (
       "intMdotVapour",
       dimensionSet(1,0,-1,0,0,0,0),
       0.0
   );

   dimensionedScalar intMdotLiquid
   (
       "intMdotLiquid",
       dimensionSet(1,0,-1,0,0,0,0),
       0.0
   );


   //cut cells with cutoff < alpha1  < 1 - cutoff and rescale
   //remaining source term
   forAll(mesh.cells(), cellI)
   {
      if(alpha1[cellI] < cutoff)
      {
          intMdotVapour.value() += alpha2[cellI]*mDot[cellI]*mesh.V()[cellI];
      }
      else if (alpha1[cellI] > 1.0 - cutoff)
      {
          intMdotLiquid.value() += alpha1[cellI]*mDot[cellI]*mesh.V()[cellI];
      }
   }

   reduce(intMdotVapour.value(),sumOp<scalar>());
   reduce(intMdotLiquid.value(),sumOp<scalar>());


   dimensionedScalar NV("NV", dimensionSet(0,0,0,0,0,0,0), 2.0);
   dimensionedScalar NL("NL", dimensionSet(0,0,0,0,0,0,0), 2.0);

   if(intMdotVapour.value() != 0){NV = evapRate/intMdotVapour;}
   if(intMdotLiquid.value() != 0){NL = evapRate/intMdotLiquid;}


   //Set source term in cells where alpha1<cutoff and alpha1>1-cutoff
   forAll(mesh.cells(), cellI)
   {
      if(alpha1[cellI] < cutoff)
      {
          continuitySource[cellI] = NV.value()*alpha2[cellI]*mDot[cellI];
      }
      else if(alpha1[cellI] > 1.0 - cutoff)
      {
          continuitySource[cellI] = - NL.value()*alpha1[cellI]*mDot[cellI];
      }
      else
      {
          continuitySource[cellI] = 0.0;
      }

   }

    Info<< "Total Evaporation Rate (before smoothing) = "
        << evapRate  << endl;
    Info<< "Total Evaporation Rate (after smoothing) = "
        << totMdot  << endl;
}
