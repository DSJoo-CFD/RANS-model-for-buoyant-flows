// Calculate Available potential energy
dimensionedScalar Tavg("Tavg",dimTemperature,gSum(mesh.V()*T.primitiveField())/gSum(mesh.V()));
Info<< "Tavg=" << Tavg.value() << "\n" << endl;

int TIndices[mesh.nCells()];
int ghIndices[mesh.nCells()];

for (int i = 0; i < mesh.nCells(); ++i)
    {
        TIndices[i] = i;
        ghIndices[i] = i;
    };

std::stable_sort(TIndices, TIndices + mesh.nCells(), [&](int a, int b) {
        return T.primitiveField()[a] < T.primitiveField()[b];});

std::stable_sort(ghIndices, ghIndices + mesh.nCells(), [&](int a, int b) {
        return gh.primitiveField()[a] > gh.primitiveField()[b];});

volScalarField TRearranged(T);
int j=0;
scalar dV = 0;
scalar dT = 0;
forAll(TRearranged, i)
{
    scalar Vi = mesh.V()[TIndices[i]];
    scalar Ti = T.internalField()[TIndices[i]];
    scalar Vj = mesh.V()[ghIndices[j]];

    while (Vj<dV) {
        TRearranged[ghIndices[j]] = dT;
        dV = dV - Vj;
        j=j+1;
        Vj = mesh.V()[ghIndices[j]];
    }

    if (Vj<dV+Vi) {
        TRearranged[ghIndices[j]] = (dV*dT + (Vj-dV)*Ti)/Vj;
        dV = Vi + dV - Vj;
        dT = Ti;
        j=j+1;
    } else {
        dT = (dT*dV + Vi*Ti)/(dV+Vi);
        dV = dV + Vi;
    }
}

// APE: Available Potential Energy per unit volume over the entire domain
dimensionedScalar APE("APE",
    dimensionSet(0,2,-2,0,0,0,0),
    beta.value()*(gSum(mesh.V()*T.primitiveField()*gh.primitiveField())
  - gSum(mesh.V()*TRearranged.primitiveField()*gh.primitiveField()))
    );

Info<< "APE=" << APE.value() << "\n" << endl;
