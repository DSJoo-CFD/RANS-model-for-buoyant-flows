Info<< "myTurbulentCorrect.H\n" << endl;

volScalarField k_(max(k,kMin_));
volScalarField epsilon_(max(epsilon,epsilonMin_));
volScalarField Tvar_(max(Tvar,TvarMin_));

volTensorField gradU(fvc::grad(U));
volVectorField gradT(fvc::grad(T));
volScalarField divU=(fvc::div(fvc::absolute(phi, U)));
volScalarField D(2.0*nu*magSqr(fvc::grad(sqrt(k))));
volScalarField D2(nu/Pr*magSqr(fvc::grad(sqrt(Tvar))));

// Damping function for k : see "LaunderSharmaKE.C" 
volScalarField fMu(exp(-3.4/sqr(scalar(1) + sqr(k_)/(nu*epsilon_)/50.0)));

// Damping function for epsilon : see "LaunderSharmaKE.C"
volScalarField f2(scalar(1) - 0.3*exp(-min(sqr(sqr(k_)/(nu*epsilon_)), scalar(50))));

// Deviatoric part of the tubulent kinetic energy production with the eddy-viscosity model : see "LaunderSharmaKE.C" 
volScalarField G=(nut*(gradU && dev(twoSymm(gradU))));


// Algebraic turbulent heat flux model
volScalarField Ceta_(min(Ceta,Ceta2*epsilon/Ct/beta/sqrt(g&g*k_*Tvar_))); // Clipping for the realizability constraint 
volScalarField stressRealizable(max(scalar(1),stessRealCoeff*0.25*nut/k_*sqrt(twoSymm(gradU)&&twoSymm(gradU))));
thf=
(
   -(Imat-Ct*k/(epsilon_+D)*Cxi*gradU)
   &
   (
       ((Imat*2/3*k-2*nut*(twoSymm(gradU)/stressRealizable))&gradT)
       +Ceta_*beta*Tvar*g
         )*Ct*k/(epsilon_+D)
);


bound(k, kMin_);

// Dissipation equation
 fvScalarMatrix epsEqn  
 (
    fvm::ddt(epsilon)
  + fvm::div(phi, epsilon)
  - fvm::laplacian(nu+nut/sigmaEps,epsilon)
 ==
    Ceps1*G*epsilon/k_
  + Cepsg*(-beta*(g&thf))*epsilon_/k_ 
  - fvm::SuSp(((2.0/3.0)*Ceps1 - Ceps3)*divU, epsilon)
  - fvm::Sp(Ceps2*f2*epsilon/k_, epsilon)
 );
epsEqn.relax();
fvOptions.constrain(epsEqn);
epsEqn.solve();
bound(epsilon, epsilonMin_);


// Turbulent kinetic energy equation
tmp<fvScalarMatrix> kEqn
(
    fvm::ddt(k)
  + fvm::div(phi, k)
  - fvm::laplacian(nu+nut/sigmaK,k)
 ==
    G - fvm::SuSp(2.0/3.0*divU, k)
  + max(-beta*(g&thf),epsilon_*scalar(0)) // Negative buoyant production has been clipped to zero
  - fvm::Sp((epsilon + D)/k_, k)
);
kEqn.ref().relax();
solve(kEqn);
bound(k, kMin_);

// Temperature variance equation
tmp<fvScalarMatrix> TvarEqn
(
   fvm::ddt(Tvar) +
   fvm::div(phi, Tvar)
 - fvm::laplacian(nu/Pr+nut/Prt,Tvar)
==
 - 2*(thf&gradT)
 - fvm::Sp(epsilon/k_/CR,Tvar)
 - D2
);
TvarEqn.ref().relax();
solve(TvarEqn);
bound(Tvar, TvarMin_);

// Update the turbulent viscosity before solving the momentum transport equation.
fMu = (exp(-3.4/sqr(scalar(1) + sqr(k)/(nu*epsilon)/50.0)));
nut = Cmu*fMu*sqr(k)/epsilon;

