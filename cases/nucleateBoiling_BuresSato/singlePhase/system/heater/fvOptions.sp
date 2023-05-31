/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2106                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvOptions;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
heatingPower
{
    type            scalarCodedSource;
    selectionMode   all;

    fields          (h);
    name            sourceTime;

    codeCorrect
        #{
            
        #};

    codeConstrain
        #{

        #};

    codeAddSup
    #{
        const scalarField& V = mesh_.V();
        scalarField& heSource = eqn.source();

        // Retrieve the x component of the cell centres
        const scalarField& cellx = mesh_.C().component(0);
	const scalar heatFlux = 481.0e03; 
	const scalar Heater_thickness = 0.5e-06;
       
        // Apply the source
        forAll(cellx, i)
        {
            // cell volume specific source
            heSource[i] -= heatFlux/Heater_thickness*(0.004-cellx[i])/0.004*V[i]; 
        };
        
    #};
}


// ************************************************************************* //
