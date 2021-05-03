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

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
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

    #include "createPhi.H"
    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );

    // Info<< "Reading field Tau\n" << endl;
    // volScalarField Tau
    // (
    //     IOobject
    //     (
    //         "Tau",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::READ_IF_PRESENT,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh,
    //     zeroGradientFvPatchScalarField::typeName
    // );
    
    // added by JX
    volSymmTensorField Tau 
    (
        IOobject
        (
            "Tau",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("Tau", dimensionSet(0,2,-2,0,0,0,0), symmTensor(0,0,0,0,0,0))
    );

    
    Tau.internalField() = 2.0/3.0*I*turbulence->k() - turbulence->nut()*twoSymm(fvc::grad(U));
    
    //volSymmTensorField Tau(
    //                        2.0/3.0*I*turbulence->k() 
    //                      - turbulence->nut()*twoSymm(fvc::grad(U))
    //                      );

    // Basic, Dimensional Field
    volTensorField UGrad(fvc::grad(U));
    volVectorField pGrad("gradP",fvc::grad(p));
    volVectorField kGrad("gradK",fvc::grad(turbulence->k()));
    volSymmTensorField S("S",symm(UGrad));
    volTensorField Omega("Omega",skew(UGrad));
    volScalarField magS(mag(S));
    volScalarField magOmega(mag(Omega));
    volScalarField magSSqr(magSqr(S));
    volScalarField magOmegaSqr(magSqr(Omega));
    volScalarField dist(wallDist(mesh).y());
    volVectorField vorticity(fvc::curl(U));
    
    volScalarField k(turbulence->k());
    volScalarField epsilon(turbulence->epsilon());

    volSymmTensorField STilda("STilda", S / (magS + epsilon / k));
    volTensorField WTilda("WTilda", Omega / (magOmega + magS));
    volVectorField pGradTilda("gradPTilda", pGrad / (mag(pGrad) + mag(U & UGrad)));
    volVectorField kGradTilda("gradKTilda", kGrad / (mag(kGrad) + epsilon / sqrt(k) * 0.01));
    

    dimensionedScalar maxReT("maxReT", dimless, 2.0);

    // Exact solution of markers
    volScalarField pGradAlongStream_Exact
    (
        IOobject
        (
            "pGradAlongStream_exact",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("pGradAlongStream_exact", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField UGradMisalignment_Exact
    (
        IOobject
        (
            "UGradMisalignment_exact",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("UGM_exact", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField QCriterion_Exact
    (
        IOobject
        (
            "QCriterion_exact",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("QC_xact", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField ReT_Exact
    (
        IOobject
        (
            "ReTurbulence_exact",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("ReT_exact", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField TurbulenceIntensity_Exact
    (
        IOobject
        (
            "TurbulenceIntensity_exact",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("TI_exact", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField TurbulenceTime_Exact
    (
        IOobject
        (
            "TurbulenceTime_exact",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("TT_exact", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField ConvectionTKE_Exact
    (
        IOobject
        (
            "ConvectionTKE_exact",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("ConvectionK_exact", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField TauRatio_Exact
    (
        IOobject
        (
            "TauRatio_exact",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("TauRatio_exact", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField PressureStress_Exact
    (
        IOobject
        (
            "PressureStress_exact",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("PressureStress_exact", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    // Validation Mode (Update the exact solution of markers)
    //#include "setSynU.H"
    //UGrad = fvc::grad(U);
    //pGrad = fvc::grad(p);
    //S = symm(UGrad);
    //Omega = skew(UGrad);
    //magSSqr = magSqr(S);
    //magOmegaSqr = magSqr(Omega);
    //magS = mag(S);
    //magOmega = mag(Omega);
    //#include "calExactSolution.H"

    // Calculation of Makers, which are nondimensionalized scalar fields

    volScalarField pGradAlongStream
    (
        IOobject
        (
            "pGradAlongStream",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("pGradAlongStream", dimensionSet(0,0,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );
    volScalarField pG_denominator(
                                    mag(U) * mag(pGrad)
                                    +
                                    mag(U & pGrad)
                                 );

    pGradAlongStream.internalField() = (U.internalField() & pGrad.internalField()) 
                                       / 
                                       Foam::max
                                       (
                                           pG_denominator,
                                           //dimensionedScalar("minpG",dimless,SMALL),
                                           dimensionedScalar("minpG",dimensionSet(0,2,-3,0,0,0,0),SMALL)
                                       );
    //pGradAlongStream.internalField() = (
    //                                       mag(U.internalField()) * mag(pGrad.internalField())
    //                                       +
    //                                       mag(U.internalField() & pGrad.internalField()) 
    //                                   );

    volScalarField UGradMisalignment
    (
        IOobject
        (
            "UGradMisalignment",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("UGradMisalignment", dimensionSet(0,0,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );
    UGradMisalignment.internalField() = mag(U.internalField() & UGrad.internalField() & U.internalField()) 
                        / 
                        ( 
                            mag(U.internalField()) * mag(UGrad.internalField() & U.internalField())
                            + 
                            mag(U.internalField() & UGrad.internalField() & U.internalField()) 
                        );

    volScalarField Lchar 
    (
        IOobject
        (
            "Lchar",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Lchar", dimensionSet(0,0,0,0,0,0,0), 0.05),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField Curvature 
    (
        IOobject
        (
            "Curvature",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Curvature", dimensionSet(0,0,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );
    Curvature.internalField() = mag(U.internalField() & UGrad.internalField())
                              / 
                              (
                                mag(U.internalField() & U.internalField())
                              + mag(U.internalField() & UGrad.internalField())
                              ); 

    volScalarField QCriterion
    (
        IOobject
        (
            "QCriterion",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("QCriterion", dimensionSet(0,0,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    QCriterion.internalField() = (magOmegaSqr.internalField() - magSSqr.internalField()) 
                                 / (magOmegaSqr.internalField() + magSSqr.internalField());

    volScalarField ReT
    (
        IOobject
        (
            "ReTurbulence",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Foam::min((Foam::sqrt(k) * dist) / (50.0 * nu), maxReT),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField TurbulenceIntensity
    (
        IOobject
        (
            "TurbulenceIntensity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        k / (0.5 * (U & U) + k),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField TurbulenceTime
    (
        IOobject
        (
            "TurbulenceTime",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("TurbulenceTime", dimensionSet(0,0,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    TurbulenceTime.internalField() = (magS.internalField() * k.internalField())
                                   / (magS.internalField() * k.internalField()
                                     + epsilon.internalField());

    volScalarField ConvectionTKE
    (
        IOobject
        (
            "ConvectionTKE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("ConvectionTKE", dimensionSet(0,0,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );
    ConvectionTKE.internalField() = (U.internalField() & kGrad.internalField())
                                    / 
                                    (
                                      mag(Tau.internalField() && S.internalField()) 
                                    + mag(U.internalField() & kGrad.internalField())
                                    );

    volScalarField TauRatio
    (
        IOobject
        (
            "TauRatio",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("TauRatio", dimensionSet(0,0,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );
    TauRatio.internalField() = mag(Tau.internalField())
                             / (k.internalField() + mag(Tau.internalField()));

    volScalarField PressureStress
    (
        IOobject
        (
            "PressureStress",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("PressureStress", dimensionSet(0,0,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );
    PressureStress.internalField() = mag(pGrad.internalField())
                                   / (
                                       mag(pGrad.internalField())
                                     + mag(3.0*cmptAv(U.internalField() & diag(UGrad.internalField())))
                                     );

    // Discretization Error of Markers
    volScalarField pGradAlongStream_Err
    (
        IOobject
        (
            "pGradAlongStream_err",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("pGradAlongStream_err", dimensionSet(0,0,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );
    //pGradAlongStream_Err.internalField() = (pGradAlongStream.internalField()
    //                                        - pGradAlongStream_Exact.internalField())
    //                                        /
    //                                        (
    //                                            mag(pGradAlongStream.internalField() 
    //                                                - pGradAlongStream_Exact.internalField())
    //                                            +
    //                                            mag(pGradAlongStream.internalField())
    //                                        );


    volScalarField UGradMisalignment_Err
    (
        IOobject
        (
            "UGradMisalignment_err",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("UGradMisalignment_err", dimensionSet(0,0,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );
    //UGradMisalignment_Err.internalField() = (UGradMisalignment.internalField()
    //                                        - UGradMisalignment_Exact.internalField())
    //                                        /
    //                                        (
    //                                            mag(UGradMisalignment.internalField() 
    //                                                - UGradMisalignment_Exact.internalField())
    //                                            +
    //                                            mag(UGradMisalignment.internalField())
    //                                        );
    volScalarField QCriterion_Err
    (
        IOobject
        (
            "QCriterion_err",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("QCriterion_err", dimensionSet(0,0,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    //QCriterion_Err.internalField() = (QCriterion.internalField()
    //                                 - QCriterion_Exact.internalField())
    //                                 /
    //                                 (
    //                                     mag(QCriterion.internalField() 
    //                                         - QCriterion_Exact.internalField())
    //                                     +
    //                                     mag(QCriterion.internalField())
    //                                 );

    // Write out the markers
    pGradAlongStream.write();
    pGradAlongStream_Exact.write();
    pGradAlongStream_Err.write();
    UGradMisalignment.write();
    UGradMisalignment_Exact.write();
    UGradMisalignment_Err.write();
    QCriterion.write();
    QCriterion_Exact.write();
    QCriterion_Err.write();
    ReT.write();
    ReT_Exact.write();
    TurbulenceIntensity.write();
    TurbulenceIntensity_Exact.write();
    TurbulenceTime.write();
    TurbulenceTime_Exact.write();
    ConvectionTKE.write();
    ConvectionTKE_Exact.write();
    TauRatio.write();
    TauRatio_Exact.write();
    PressureStress.write();
    PressureStress_Exact.write();
    Curvature.write();
    // Added by JX
    Tau.write();
    

    //output tensor and vector bases for the generation of more features
    pGradTilda.write();
    kGradTilda.write();
    STilda.write();
    WTilda.write();

    // Initialize region flag
    volScalarField regionFlag
    (
        IOobject
        (
            "regionFlag",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("regionFlag", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    // Set region flag
    #include "regionInfo.H"
