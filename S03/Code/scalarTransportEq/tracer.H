fvScalarMatrix TracerEqn
(
	fvm::ddt(Tracer) + fvm::div(phi, Tracer) - fvm::laplacian(DT, Tracer)
//	==
//	fvModels.source(T)
);

TracerEqn.solve();

Info<< "Min, Max Tracer:" <<  min(Tracer).value() <<  ", " <<  max(Tracer).value() << endl;
