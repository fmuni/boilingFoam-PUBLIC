Info << "Updating meshes" << endl;
for(int meshIt = 0; meshIt < nMeshCorr; meshIt++)
{
    //- Update mesh in solid and fluid to synchronize AMR
    forAll(fluidRegions,i)
    {
        fluidRegions[i].update();

        //- Correct interface upon refinement        
        if(fluidRegions[i].changing())
        {
            interfaceFluid[i].correct();
        }
    }

    forAll(solidRegions,i)
    {
        //- Fetch field from nbr mesh
        refineSolid[i].correctBoundaryConditions();
     
        solidRegions[i].update();
    }
}
