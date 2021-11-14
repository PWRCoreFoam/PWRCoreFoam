
    Foam::Info
        << "Create subChannel mesh for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

    Foam::subChannelMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    mesh.buildSubChannelMesh();


    mesh.output_mesh_information();

    mesh.createRods();

    mesh.modify_mesh();
    mesh.modify_mesh_for_core_edge();
    mesh.modify_mesh_for_assembly_edge();

    mesh.output_mesh_information();

    surfaceScalarField * magSfPtr_ = mesh.get_magSfPtr();
    surfaceScalarField & faces_area_mag = *magSfPtr_;


