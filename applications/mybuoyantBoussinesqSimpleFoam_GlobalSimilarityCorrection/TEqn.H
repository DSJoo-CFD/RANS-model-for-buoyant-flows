Info<< "TEqn.H\n" << endl;


D=(2.0*nu*magSqr(fvc::grad(sqrt(k))));
k_ = max(k,kMin_);
epsilon_ = max(epsilon+D,epsilonMin_);
volScalarField alphaEff("alphaEff", nu/Pr + Ct*sqr(k)/epsilon_);

Ceta_ = (min(Ceta,Ceta2*epsilon/Ct/beta/sqrt(g&g*k_*Tvar_)));
stressRealizable = (max(scalar(1),stessRealCoeff*0.25*nut/k_*sqrt(twoSymm(gradU)&&twoSymm(gradU))));

// thf2: The subtraction of the diffusion-like part from the algebraic turbulent heat flux model. 
//    The algebraic heat flux model has been decomposed into the diffusion-like part 
//    (which is equal to the gradient diffusion model) and the rest of the terms. 
//    Then, the diffusion-like part is solved with the Laplacian scheme (which is designed to solve diffusion) in OpenFOAM.
//    The inverse matrix part of the algebraic model, (I+T)^{-1}, has been approximated as (I-T) by Neumann series.
thf2=
(
   - Ct*k/epsilon_*(Ceta_*beta*Tvar*g-nut*(twoSymm(gradU)&gradT)/stressRealizable)
   + pow(Ct*k/epsilon_,2)*Cxi*(gradU&
      (
         ((Imat*2/3*k-nut*(twoSymm(gradU)))&gradT/stressRealizable)
          +Ceta_*beta*Tvar*g
      ))
);

// Eliminate checkerboard oscillations (odd-even decoupling) in the temperature field using face interpolation
surfaceVectorField thf2f("thf2f", fvc::interpolate(thf2));
thf2 = fvc::average(thf2f); 

fvScalarMatrix TEqn
(
    fvm::div(phi, T)
  - fvm::laplacian(alphaEff, T)
  + fvc::div(thf2)
 ==
    radiation->ST(rhoCpRef, T)
  + Q_
  + fvOptions(T)
);

TEqn.relax();

fvOptions.constrain(TEqn);

TEqn.solve();

radiation->correct();

fvOptions.correct(T);

rhok = 1.0 - beta*(T - TRef);
