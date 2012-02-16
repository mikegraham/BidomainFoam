// Copyright 2012 Michael Graham

#include "fitzHughNagumo.H"

#include <math.h>

Foam::fitzHughNagumo::fitzHughNagumo(
    scalar epsilon,
    scalar gamma,
    scalar beta
)
{
    this->epsilon_ = epsilon;
    this->gamma_ = gamma;
    this->beta_ = beta;
}

Foam::label Foam::fitzHughNagumo::nEqns() const
{
    return 2;
}

void Foam::fitzHughNagumo::derivatives
(
    const Foam::scalar t,
    const Foam::scalarField& y,
    Foam::scalarField& dydt
) const
{
    dydt[0] = (1 / this->epsilon_) * (y[0] - pow(y[0], 3) - y[1]);
    dydt[1] = this->epsilon_ * (y[0] - this->gamma_ * y[1] + this->beta_);
}

void Foam::fitzHughNagumo::jacobian
(
    const Foam::scalar t,
    const Foam::scalarField& y,
    Foam::scalarField& dfdt,
    Foam::scalarSquareMatrix& dfdy
) const
{
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;

    dfdy[0][0] = (1 - pow(y[0], 2)) / this->epsilon_;
    dfdy[0][1] = -1.0 / this->epsilon_;
    dfdy[1][0] = this->epsilon_ * y[0];
    dfdy[1][1] = -this->epsilon_ * this->gamma_;
}

