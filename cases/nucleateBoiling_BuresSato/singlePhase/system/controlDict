/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  plus                                  |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     icoBoilingFoam;

startFrom       latestTime;

startTime       0;

stopAt          endTime;

endTime         0.0881;

deltaT         1e-10;

writeControl    adjustableRunTime;

writeInterval   1e-3;

purgeWrite      2;

writeFormat     binary;

writePrecision  10;

writeCompression no;

timeFormat      general;

timePrecision   10;

runTimeModifiable yes;

adjustTimeStep  yes;

maxCo           0.25;
maxDi           10000;
maxAlphaCo      0.25;
maxDeltaT       1.0;//#calc "1e-3*$t0";;
maxCapillaryNum 0.5;

functions
{

//----------------------------------------------------------------------------//
//- Function to sample the interface.
        freeSurface
   {
       type            surfaces;
       libs            ("libsampling.so");
       surfaceFormat  raw;
       executeControl  timeStep;
       executeInterval 500;
       writeControl    timeStep;
       writeInterval 500;
       sampleOnExecute true;
       region fluid;
       fields
       (
           alpha.liquid
       );
       surfaces
       (
           freeSurface
           {
               type        isoSurfaceCell;
               isoField    alpha.liquid;
               isoValue    0.5;
               interpolate true;
               regularise  false;
           }
       );
       interpolationScheme cellPointFace;
   }

bubbleVol
    {
        type            volFieldValue;
        libs            (fieldFunctionObjects);
        writeControl    timeStep;
        writeInterval   10;
        writeFields     false;
        log             true;
	region 		fluid;
        operation       volIntegrate;
        fields 		(alpha.vapour);
    }

sampled_heater
       {
          type    surfaces;
          libs    (sampling);
          log     false;

          executeControl  runTime;
          executeInterval 0.001;
          writeControl    runTime;
          writeInterval 0.001;
          sampleOnExecute true;

          surfaceFormat   raw;
          region heater;
          fields      (T );

           surfaces
           {
                   heater
                   {
                           type patch;
                           patches (heater_to_fluid);
                   }
           }
         }


}
// ************************************************************************* //
