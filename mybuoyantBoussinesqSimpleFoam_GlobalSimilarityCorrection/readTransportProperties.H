singlePhaseTransportModel laminarTransport(U, phi);

// Thermal expansion coefficient [1/K]
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    laminarTransport
);

// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

// Laminar Prandtl number
dimensionedScalar Pr("Pr", dimless, laminarTransport);

// Kinematic viscosity
dimensionedScalar nu("nu",dimViscosity,laminarTransport);

// Heat source per unit volume and unit heat capacity [K/s]
scalar Q(laminarTransport.getOrDefault<scalar>("Q", 0)); 
dimensionedScalar Q_(dimensionSet(0,0,-1,1,0,0,0),Q);


