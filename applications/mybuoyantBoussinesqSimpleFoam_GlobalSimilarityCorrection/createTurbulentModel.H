/*---------------------------------------------------------------------------*\
Models for k, epsilon, and nut are implemented similarly to those in 
"LaunderSharmaKE.C" in OpenFOAM. The algebraic turbulent heat flux model and 
the temperature variance equation are implemented as described in the paper by 
Kenjereš, S. and Hanjalić, K., 1995, Int. J. Heat Fluid Flow, 16(5).

My new modeling approach only algebraically modifies the buoyancy-related term 
in the epsilon equation by incorporating a global effect of buoyancy, 
represented by available potential energy.
\*---------------------------------------------------------------------------*/


// Turbulent viscosity model
dimensionedScalar Cmu(laminarTransport.getOrDefault<scalar>("Cmu", 0.09));


// Model transport equation for k (turbulent kinetic energy)
Info<< "Reading field k\n" << endl;
volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
dimensionedScalar sigmaK(laminarTransport.getOrDefault<scalar>("sigmaK", 1.0));


// Model transport equation for epsilon (dissipation rate)
Info<< "Reading field epsilon\n" << endl;
volScalarField epsilon
(
    IOobject
    (
        "epsilon",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
dimensionedScalar Ceps1(laminarTransport.getOrDefault<scalar>("Ceps1", 1.44));
dimensionedScalar Ceps2(laminarTransport.getOrDefault<scalar>("Ceps2", 1.92));
dimensionedScalar Ceps3(laminarTransport.getOrDefault<scalar>("Ceps3", -0.33));
dimensionedScalar sigmaEps(laminarTransport.getOrDefault<scalar>("sigmaEps", 1.3));


// Turbulent viscosity
volScalarField nut
(
    IOobject
    (
        "nut",
        runTime.timeName(),
        mesh
    ),
    Cmu*sqr(k)/epsilon
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Algebriac turbulent heat flux model
volTensorField Imat
(
    IOobject
    (
        "Imat",
        runTime.timeName(),
        mesh
    ),
    mesh,
    tensor(1,0,0,0,1,0,0,0,1)
);
dimensionedScalar Ct(laminarTransport.getOrDefault<scalar>("Ct", 0.15));
dimensionedScalar Cxi(laminarTransport.getOrDefault<scalar>("Cxi", 0.6));
dimensionedScalar Ceta(laminarTransport.getOrDefault<scalar>("Ceta", 0.6));
dimensionedScalar Ceta2(laminarTransport.getOrDefault<scalar>("Ceta2", 0.2)); 
    //Ceta2: A coefficient associated with the clipping treatment to ensure realizability of the algebraic turbulent heat flux model
dimensionedScalar stessRealCoeff(laminarTransport.getOrDefault<scalar>("stessRealCoeff", 5)); 
    //stressRealCoeff: Clip Reynolds stress in the algebraic turbulent heat flux model to satisfy the realizability condition. 
    //      Set an empirical value equal to or larger than one. 
    //      A larger value makes the model more stable, but there are a few risks of damping Reynolds stress unphysically.



// New modeling coefficients modify the buoyancy-related term in the epsilon equation.
dimensionedScalar Cepsg(laminarTransport.getOrDefault<scalar>("Cepsg", 0.8));
dimensionedScalar Ca(laminarTransport.getOrDefault<scalar>("Ca", 0.75));
dimensionedScalar Cb(laminarTransport.getOrDefault<scalar>("Cb", 0.5));
dimensionedScalar Cd(laminarTransport.getOrDefault<scalar>("Cd", 0.5));
dimensionedScalar Ce(laminarTransport.getOrDefault<scalar>("Ce", 0.4));


// Model transport equation for Tvar (temperature variance; i.e., Reynolds-averaged value of squared temperature fluctuations)
Info<< "Reading field Tvar(temperature variance)\n" << endl;
volScalarField Tvar
(
    IOobject
    (
        "Tvar",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
dimensionedScalar Prt(laminarTransport.getOrDefault<scalar>("Prt", 1.0));
dimensionedScalar CR(laminarTransport.getOrDefault<scalar>("CR", 0.75)); 
    // CR: ratio of mechanical to thermal turbulent time scales, used in the dissipation rate of temperature variance.
    // The dissipation rate of temperature variance is algebraically modeled as "(1/CR) * (epsilon / k) * Tvar".


scalar kMin(laminarTransport.getOrDefault<scalar>("kMin",1.e-13));
scalar epsilonMin(laminarTransport.getOrDefault<scalar>("epsilonMin",1.e-10));
scalar TvarMin(laminarTransport.getOrDefault<scalar>("TvarMin", 1.e-13));
dimensionedScalar kMin_(dimensionSet(0,2,-2,0,0,0,0),kMin);
dimensionedScalar epsilonMin_(dimensionSet(0,2,-3,0,0,0,0),epsilonMin);
dimensionedScalar TvarMin_(dimensionSet(0,0,0,2,0,0,0),TvarMin);



// Turbulent heat flux
volVectorField thf
(
    IOobject
    (
        "thf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -Ct*Ceta*g*beta*k/epsilon*Tvar
);


volVectorField Nu
(
    IOobject
    (
        "Nu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -Ct*Ceta*g*beta*k/epsilon*Tvar
);

volVectorField thf2
(
    IOobject
    (
        "thf2",
        runTime.timeName(),
        mesh
    ),
    -Ct*Ceta*g*beta*k/epsilon*Tvar
);


Info<< "Transport properties:" << endl;
Info<< "nu = " << nu.value() << endl;
Info<< "beta = " << beta.value() << endl;
Info<< "Q = " << Q << endl;
Info<< "Pr = " << Pr.value() << endl;
dimensionedScalar Rai(mag(g)*beta*Q/nu/nu/nu*Pr*Pr);
Info<< "Rai = " << Rai.value() << endl;
Info<< " " << endl;


Info<< "Modeling coefficients:" << endl;
Info<< "Ceps1 = " << Ceps1.value() << endl;
Info<< "Ceps2 = " << Ceps2.value() << endl;
Info<< "Ceps3 = " << Ceps3.value() << endl;
Info<< "Cmu = " << Cmu.value() << endl;

Info<< "Ct = " << Ct.value() << endl;
Info<< "Cxi = " << Cxi.value() << endl;
Info<< "Ceta = " << Ceta.value() << endl;
Info<< "CR = " << CR.value() << endl;

Info<< "sigmaK = " << sigmaK.value() << endl;
Info<< "sigmaEps = " << sigmaEps.value() << endl;
Info<< "Prt = " << Prt.value() << endl;

Info<< "Cepsg = " << Cepsg.value() << endl;
Info<< "Ca = " << Ca.value() << endl;
Info<< "Cb = " << Cb.value() << endl;
Info<< "Cd = " << Cd.value() << endl;
Info<< "Ce = " << Ce.value() << "\n" << endl;


