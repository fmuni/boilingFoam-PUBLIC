//- Set alpha between 0 and 1 only at the interface. Remove "dust".
{

    Switch hardBoundAlpha
    (
        mixture.lookupOrDefault<Switch>("alphaInterfaceOnly",false)
    );

    if(hardBoundAlpha)
    {
        const boolList& interfaceCells(interface.surf().interfaceCell());

        forAll(interfaceCells,celli)
        {
            if(!interfaceCells[celli])
            {
                alpha1[celli] = pos0(alpha1[celli]-0.5);
            }
        }
    }

}
