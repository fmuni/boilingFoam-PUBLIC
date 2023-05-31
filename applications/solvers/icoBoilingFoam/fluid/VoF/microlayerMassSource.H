//- Update mDot to account for microlayer models at the walls

if(mixture.found("microlayerModel"))
{
    //Info<< "\nStarting microlayer stuff"<< endl;
    dictionary microDict(mixture.subDict("microlayerModel"));
    Switch debugML(mixture.lookupOrDefault<Switch>("debug",false));

    //- Coefficients for power law
    // dimensionedScalar b
    // (
    //     "b",
    //     dimless,
    //     microDict
    // );
    //
    // //-Note that the dimensions of A are based on a reference temperature =1K.
    // // This to avoid awkward dimensions due to the exponential of T-Tsat.
    // dimensionedScalar A
    // (
    //     "A",
    //     dimPower/dimLength,
    //     microDict
    // );

    //- Coefficients for polynomium
    scalarField coeffs(microDict.lookup("coefficients"));

    if(coeffs.size() != 4)
    {
        Info<<"\nThe microlayer model requires 4 coefficients!\n"
            << abort(FatalError);
    }

    wordList activePatches(microDict.lookup("patches"));

    scalarField ap(volPointInterpolation::New(mesh).interpolate(alpha1));

    cutFaceIso faceCut(mesh,ap);

    forAll(activePatches,patchi)
    {
        label patchID(mesh.boundaryMesh().findPatchID(activePatches[patchi]));

        if(patchID == -1)
        {
            continue;
        }

        const fvPatchScalarField& TBf(T.boundaryFieldRef()[patchID]);
        const fvPatchScalarField& aBf(alpha1.boundaryFieldRef()[patchID]);
        scalarField& kappaBf(kappa.boundaryFieldRef()[patchID]);
        const labelUList& fCells(TBf.patch().faceCells());
        const scalarField& deltas(TBf.patch().deltaCoeffs());
        const scalarField pA(mesh.magSf().boundaryField()[patchID]);

        //- Cache wall gradient
        scalarField dTdn(TBf.snGrad());

        const label startFace(TBf.patch().start());

        forAll(aBf,facei)
        {

            label cellI = fCells[facei];

            mDot0[cellI] = 0;
            energySource[cellI] = 0;
            energyExp[cellI] = 0;


            //- Select interface cells only
            if(aBf[facei] < 1e-3 || aBf[facei] > 0.999)
            {
                continue;
            }

            if(faceCut.calcSubFace(startFace + facei,0.5) != 0)
            {
                continue;
            }

            // if(faceCut.subFacePoints().size() == 0)
            // {
            //     continue;
            // }
            //Info<< "\npoints = "<< faceCut.surfacePoints() <<  endl;

            //- Compute the contact length
            point p0(faceCut.surfacePoints()[0]);
            point p1(faceCut.surfacePoints()[1]);

            dimensionedScalar lc
            (
                "lc",
                dimLength,
                mag(p0-p1)
            );

            //- Set conduction to zero.
            //  This does not affect cht, which is based on the true kappa.
            kappaBf[facei] = 0.;

            scalar deltaTmp(TBf[facei] - Tsat.value());

            //- Re-compute mDot0 at the interface
            mDot0[cellI] =
            (
                lc.value()/(hf.value()*mesh.V()[cellI])
              * (
                    (coeffs[0]*deltaTmp)
                  + (coeffs[1]*pow(deltaTmp,2))
                  + (coeffs[2]*pow(deltaTmp,3))
                  + (coeffs[3]*pow(deltaTmp,4))
                )
            );

            energySource[cellI] =
            (
                0
            );

            energyExp[cellI] =
            (
                -lc.value()/mesh.V()[cellI]
              * (
                    (coeffs[0]*(deltaTmp))
                  + (coeffs[1]*pow(deltaTmp,2))
                  + (coeffs[2]*pow(deltaTmp,3))
                  + (coeffs[3]*pow(deltaTmp,4))
                )
            );

            if(debugML)
            {
                fvPatchScalarField& lcBf
                (
                    lcFluid[i].boundaryFieldRef()[patchID]
                );

                lcBf[facei] = lc.value();
            }

            // energySource[cellI] =
            // (
            //     -    lc.value()*coeffs[0]/mesh.V()[cellI]
            // );
            //
            // energyExp[cellI] =
            // (
            //     -lc.value()/mesh.V()[cellI]
            //   * (
            //         (coeffs[0]*(Tsat.value() - deltas[facei]*dTdn[facei]))
            //       + (coeffs[1]*pow(deltaTmp,2))
            //       + (coeffs[2]*pow(deltaTmp,3))
            //       + (coeffs[3]*pow(deltaTmp,4))
            //     )
            // );

            //- Compute coefficient
            // dimensionedScalar coeffA
            // (
            //     A*lc
            // );
            //
            // //- No layer if the temperature is close to Tsat
            // if(TBf[facei] - Tsat.value() < SMALL)
            // {
            //     continue;
            // }
            //
            // //- Re-compute mDot0 at the interface
            // mDot0[cellI] =
            // (
            //     coeffA.value()
            //   * Foam::pow
            //     (
            //         TBf[facei] - Tsat.value(),
            //         b.value()
            //     )
            // );
            //
            // //- Compute the linearized coefficient of the energy equation
            // //  assuming that T_p = T_c - dT/dn|_p (xc-x_p)
            //
            // energySource[cellI] =
            // (
            //     coeffA.value()
            //   * Foam::pow
            //     (
            //         TBf[facei] - Tsat.value(),
            //         b.value() - 1.
            //     )
            // );
            //
            // energyExp[cellI] =
            // (
            //     energySource[cellI]
            //   * (
            //         Tsat.value() - deltas[facei]*dTdn[facei]
            //     )
            // );
        }
    }



}
