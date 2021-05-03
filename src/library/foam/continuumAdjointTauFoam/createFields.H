    Info << "Reading transportProperties\n" << endl;
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

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

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );

    dimensionedScalar unitTime("unitTime", dimTime, 1.0); 
    dimensionedScalar unitLength("unitLength", dimLength, 1.0); 
    dimensionedScalar smallNut("smallNut", dimLength*dimLength/dimTime, 1e-10); 

    Info<< "Reading field U\n" << endl;

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField UDNS
    (
        IOobject
        (
            "UDNS",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );



    #include "createPhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);


    Info<< "Reading field pa\n" << endl;
    volScalarField pa
    (
        IOobject
        (
            "pa",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field Ua\n" << endl;
    volVectorField Ua
    (
        IOobject
        (
            "Ua",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhia.H"


    label paRefCell = 0;
    scalar paRefValue = 0.0;
    setRefCell
    (
        pa,
        mesh.solutionDict().subDict("SIMPLE"),
        paRefCell,
        paRefValue
    );


    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );


    dimensionedTensor zeroSensitivity("0", dimVelocity/dimLength, tensor::zero);
    dimensionedTensor zeroTau("0", dimVelocity*dimVelocity, tensor::zero);
    dimensionedScalar zeroNut("0", dimVelocity*dimLength, scalar(0.0));

    dimensionedScalar lambda(laminarTransport.lookup("lambda"));
    dimensionedScalar regulationCoeff(laminarTransport.lookup("regulationCoeff"));
    dimensionedVector explicitGradient(laminarTransport.lookup("explicitGradient"));

    label nSimpleLoop = laminarTransport.lookupOrDefault("nSimpleLoop",1);
    Info << "nSimpleLoop is: " << nSimpleLoop << endl;

    label nAdjointLoop = laminarTransport.lookupOrDefault("nAdjointLoop",1);
    Info << "nAdjointLoop is: " << nAdjointLoop << endl;

    word adjointType = laminarTransport.lookup("adjointType");

    // const labelList& inletCells = mesh.boundary()["inlet"].faceCells();
    // const labelList& outletCells = mesh.boundary()["outlet"].faceCells();

    volTensorField tauR
    (
        IOobject
        (
            "tauR",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        zeroTau
    );
    //zeroCells(Tau, inletCells);
    //zeroCells(Tau, outletCells);
    
    volScalarField nut
    (
        IOobject
        (
            "nut",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        zeroNut
    );

    volScalarField l2nNut("l2nNut", fvc::grad(nut) & fvc::grad(nut));

    volTensorField sensitiveJ
    (
        IOobject
        (
            "sensitiveJ",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        zeroSensitivity
    );

    dimensionedScalar maximumJ(laminarTransport.lookup("maximumJ"));
