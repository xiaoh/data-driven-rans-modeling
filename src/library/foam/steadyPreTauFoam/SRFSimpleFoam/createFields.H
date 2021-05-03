    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field Urel\n" << endl;
    volVectorField Urel
    (
        IOobject
        (
            "Urel",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading/calculating face flux field phi\n" << endl;
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(Urel) & mesh.Sf()
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);

    singlePhaseTransportModel laminarTransport(Urel, phi);

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(Urel, phi, laminarTransport)
    );

    Info<< "Creating SRF model\n" << endl;
    autoPtr<SRF::SRFModel> SRF(SRF::SRFModel::New(Urel));
