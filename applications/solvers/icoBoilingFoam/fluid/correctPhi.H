CorrectPhi
(
    U,
    phi,
    p_rgh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU)),
    absDivU,
    pimple
);

#include "continuityErrs.H"
