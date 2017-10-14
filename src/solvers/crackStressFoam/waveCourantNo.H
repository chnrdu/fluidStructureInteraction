{
    scalar E0 = stress().rheology().law().E()().internalField()[0];
    scalar nu0 = stress().rheology().law().nu()().internalField()[0];
    scalar rho0 = stress().rheology().rho()().internalField()[0];

    scalar waveVelocity =
        ::sqrt(E0*(1 - nu0)/(rho0*(1 + nu0)*(1 - 2*nu0)));

    // Courant number
    scalarField Co =
        waveVelocity*runTime.deltaT().value()
       *mesh.surfaceInterpolation::deltaCoeffs().internalField();

    Info<< "Courant Number mean: " << gAverage(Co)
        << " max: " << gMax(Co)
        << " wave velocity magnitude: " << waveVelocity << endl;
}
