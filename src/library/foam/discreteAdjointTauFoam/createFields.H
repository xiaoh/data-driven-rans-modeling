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

Info<< "Reading field UDNS\n" << endl;
volVectorField UDNS
(
    IOobject
    (
        "UDNS",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);
autoPtr<incompressible::turbulenceModel> turbulence
(
	incompressible::turbulenceModel::New(U, phi, laminarTransport)
);
// Towara: create field to hold resistance term tau
dimensionedSymmTensor zeroSensitivity("0", dimVelocity*dimVelocity, symmTensor::zero);
dimensionedSymmTensor zeroTau("0", dimVelocity*dimVelocity, symmTensor::zero);

volSymmTensorField tauR
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
    zeroTau,
    zeroGradientFvPatchScalarField::typeName
);

volSymmTensorField sens
(
    IOobject
    (
        "sens",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    zeroSensitivity,
    zeroGradientFvPatchScalarField::typeName
);
