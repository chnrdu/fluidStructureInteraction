Info << " Reading dynamic fluid2 mesh" << endl;

fvMesh fluidMesh2
(
    IOobject
    (
        "fluid2",
//         runTime.constant(),
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
    )
);

