// Copyright 2012 Michael Graham

#include "fvCFD.H"
#include "ODE.H"
#include "ODESolver.H"
#include "RK.H"
#include "fitzHughNagumo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    
    const Foam::fitzHughNagumo ode(0.2, 0.8, 0.7);
    Foam::autoPtr<Foam::ODESolver> odeSolver = Foam::ODESolver::New("RK", ode);
    scalarField odeY(ode.nEqns());
    odeY[0] = 0.0;
    odeY[1] = 0.0;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // Solve the ODEs
        scalar tEnd = runTime.timeOutputValue();
        scalar dT = runTime.deltaT().value();
        scalar tStart = tEnd - dT;
        scalar eps(0.00001);
        
        Info << "Solving ODE from " << tStart << "to" << tEnd << endl;
        Info << odeY << endl;
        odeSolver->solve(
            ode, 
            tStart, 
            tEnd, 
            odeY, 
            eps, 
            dT
        );
        Info << odeY << endl;
        
        
        // FIXME: some constants sitting in for real values
        dimensionSet AmperesPerCubicMeter(0, -3, 0, 0, 0, 1, 0);
        dimensionSet AmperesPerSquareMeter(0, -2, 0, 0, 0, 1, 0);
        dimensionSet siemensPerMeter(-1, -3, 3, 0, 0, 2, 0);

        dimensionedScalar Iion(
            "Iion",
            AmperesPerSquareMeter,
            1.0
        );
        dimensionedScalar Ie(
            "Ie", 
            AmperesPerCubicMeter,
            1.0
        ); 
        dimensionedScalar sigmaBarI(
            "sigmaBarI", 
            siemensPerMeter,
            1.0
        );
        dimensionedScalar sigmaBarSum(
            "sigmaBarSum", 
            siemensPerMeter,
            1.0
        );      
        
        
        // Solve the parabolic PDE.      
        solve(
            betam * Cm * fvm::ddt(Vm)
          - fvm::laplacian(sigmai, Vm)
          - fvc::laplacian(sigmai, phie)
          + betam * Iion
        );

        // Solve the elliptic PDE
        solve
        (
            fvm::laplacian(sigmaBarSum, phie)
          + fvc::laplacian(sigmai, Vm)
          + Ie
        );

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
