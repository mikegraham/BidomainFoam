// Copyright 2012 Michael Graham

#ifndef fitzHughNagumo_H
#define fitzHughNagumo_H

#include "ODE.H"

namespace Foam {

class fitzHughNagumo
:
    public ODE
{
public:
    // Basic constructor
    fitzHughNagumo(scalar epsilon, scalar gamma, scalar beta);
    
    // Override the members defined by ODE.
    label nEqns() const;

    void derivatives
    (
        const scalar time,
        const scalarField& y,
        scalarField& dydt
    ) const;

    void jacobian
    (
        const scalar t,
        const scalarField& y,
        scalarField& dfdt,
        scalarSquareMatrix& dfdy
    ) const;
    
private:
    scalar epsilon_;
    scalar gamma_;
    scalar beta_;
};

}   

#endif
