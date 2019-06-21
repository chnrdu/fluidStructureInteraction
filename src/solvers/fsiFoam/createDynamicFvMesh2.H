Info<< "Create dynamic mesh2 for time = "
        << runTime.timeName() << nl << endl;

    autoPtr<dynamicFvMesh> mesh2Ptr
    (
        dynamicFvMesh::New
        (
            IOobject
           (

               "fluid2",
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );

    dynamicFvMesh& mesh2 = mesh2Ptr();


// // Correct symmetryPlane points

// vectorField& newPoints = const_cast<vectorField&>(stressMesh.points());

// forAll(stressMesh.boundaryMesh(), patchI)
// {
//     if (isA<symmetryPolyPatch>(stressMesh.boundaryMesh()[patchI]))
//     {
//         const labelList& meshPoints =
// 	    stressMesh.boundaryMesh()[patchI].meshPoints();

// 	vector avgN =
// 	    gAverage(stressMesh.boundaryMesh()[patchI].pointNormals());

// 	vector i(1, 0, 0);
// 	vector j(0, 1, 0);
// 	vector k(0, 0, 1);

// 	if (mag(avgN&i) > 0.95)
// 	{
// 	    forAll(meshPoints, pI)
// 	    {
// 	        newPoints[meshPoints[pI]].x() = 0;
// 	    }
// 	}
// 	else if (mag(avgN&j) > 0.95)
// 	{
// 	    forAll(meshPoints, pI)
// 	    {
// 	        newPoints[meshPoints[pI]].y() = 0;
// 	    }
// 	}
// 	else if (mag(avgN&k) > 0.95)
// 	{
// 	    forAll(meshPoints, pI)
// 	    {
// 	        newPoints[meshPoints[pI]].z() = 0;
// 	    }
// 	}

// 	avgN = gAverage(stressMesh.boundaryMesh()[patchI].pointNormals());

// 	Pout << "avgN = " << avgN << endl;
//     }
// }

// stressMesh.movePoints(newPoints);
// stressMesh.V00();
// stressMesh.moving(false);

