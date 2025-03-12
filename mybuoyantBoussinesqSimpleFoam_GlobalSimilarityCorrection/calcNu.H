// For convenience, the reference length and temperature difference scales are set to one.
{
    Nu = (T*U+thf-nu/Pr*gradT);
    dimensionedVector Nuavg("Nuavg",dimless,gSum(mesh.V()*Nu.primitiveField())/gSum(mesh.V())/nu.value()*Pr.value());

    Info<< "Nu_avg=" << Nuavg.value() << "\n" << endl;
}

