{
    forAll (stress->solution().boundaryField(), patchI)
    {
        if
        (
            stress->solution().boundaryField()[patchI].type() 
         == cohesiveZoneIncrementalFvPatchVectorField::typeName
        )
        {
            cohesiveZoneIncrementalFvPatchVectorField& cohesivePatch = 
                refCast<cohesiveZoneIncrementalFvPatchVectorField>
                (
                    stress->solution().boundaryField()[patchI]
                );

            nCrackedFaces = cohesivePatch.updateCrack();	    

            label nFaces = 
                gSum(cohesivePatch.crackIndicator())
              + gSum(cohesivePatch.crazeIndicator());

            Info << "nCrackedFaces = " << nFaces << endl;
        }
    }

    reduce(nCrackedFaces, sumOp<label>());
}
