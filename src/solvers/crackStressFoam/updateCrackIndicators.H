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

            crack.boundaryField()[patchI] = cohesivePatch.crackIndicator();
            craze.boundaryField()[patchI] = cohesivePatch.crazeIndicator();
            relSepDist.boundaryField()[patchI] = 
                cohesivePatch.relSeparationDistance();
        }
    }
}
