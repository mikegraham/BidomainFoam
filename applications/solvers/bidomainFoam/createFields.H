
    Info<< "Reading bidomainProperties" << endl;

    IOdictionary bidomainProperties
    (
        IOobject
        (
            "bidomainProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar Cm
    (
        bidomainProperties.lookup("Cm")
    );
    
    dimensionedScalar betam
    (
        bidomainProperties.lookup("betam")
    );
    
    dimensionedTensor sigmai
    (
        bidomainProperties.lookup("sigmai")
    );

    Info<< "Reading field Vm\n" << endl;
    volScalarField Vm
    (
        IOobject
        (
            "Vm",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field phie\n" << endl;
    volScalarField phie
    (
        IOobject
        (
            "phie",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
