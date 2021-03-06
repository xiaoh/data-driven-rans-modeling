/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::incompressible::RASModels::myLRR

Group
    grpIcoRASTurbulence

Description
    Launder, Reece and Rodi Reynolds-stress turbulence model for
    incompressible flows.

    The default model coefficients correspond to the following:
    \verbatim
        myLRRCoeffs
        {
            Cmu         0.09;
            Clrr1       1.8;
            Clrr2       0.6;
            C1          1.44;
            C2          1.92;
            Cs          0.25;
            Ceps        0.15;
            sigmaEps    1.3;
            couplingFactor  0.0;    // only for incompressible
        }
    \endverbatim

SourceFiles
    myLRR.C

\*---------------------------------------------------------------------------*/

#ifndef myLRR_H
#define myLRR_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class myLRR Declaration
\*---------------------------------------------------------------------------*/

class myLRR
:
    public RASModel
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;

            dimensionedScalar Clrr1_;
            dimensionedScalar Clrr2_;

            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar Cs_;
            dimensionedScalar Ceps_;
            dimensionedScalar sigmaEps_;

            dimensionedScalar couplingFactor_;


        // Fields

            //volSymmTensorField R_;
            volSymmTensorField Tau_;
            volScalarField k_;
            volScalarField epsilon_;
            volScalarField nut_;


public:

    //- Runtime type information
    TypeName("myLRR");

    // Constructors

        //- Construct from components
        myLRR
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~myLRR()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for R
        tmp<volScalarField> DREff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DREff", nut_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_/sigmaEps_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const
        {
            return Tau_;
        }

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
