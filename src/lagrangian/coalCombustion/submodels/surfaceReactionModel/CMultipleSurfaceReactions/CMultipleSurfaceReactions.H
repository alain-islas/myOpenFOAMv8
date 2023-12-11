/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::CMultipleSurfaceReactions

Description
    Kinetic/diffusion limited rate surface reaction model for coal parcels.
    Aplicable to:

    R1 C(s) + O2 -> CO2
    R2 C(s) + CO2 -> 2CO
    R3 C(s) + H2O -> H2 + CO

\*---------------------------------------------------------------------------*/

#ifndef CMultipleSurfaceReactions_H
#define CMultipleSurfaceReactions_H

#include "SurfaceReactionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
template<class CloudType>
class CMultipleSurfaceReactions;

/*---------------------------------------------------------------------------*\
            Class CMultipleSurfaceReactions Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CMultipleSurfaceReactions
:
    public SurfaceReactionModel<CloudType>
{
    // Private Data

        // Model constants

            //- Mass diffusion limited rate constant, R1 C1
            const scalar C11_;
            //- Kinetics limited rate pre-exponential constant, R1 C2
            const scalar C21_;
            //- Kinetics limited rate activation energy, R1
            const scalar E1_;

           //- Mass diffusion limited rate constant, R2 C1
            const scalar C12_;
            //- Kinetics limited rate pre-exponential constant, R2 C2
            const scalar C22_;
            //- Kinetics limited rate activation energy, R2
            const scalar E2_;

           //- Mass diffusion limited rate constant, R3 C1
            const scalar C13_;
            //- Kinetics limited rate pre-exponential constant, R3 C2
            const scalar C23_;
            //- Kinetics limited rate activation energy, R3
            const scalar E3_;


        // Addressing

            //- Cs positions in global/local lists
            label CsLocalId_;
            //- O2 position in global list
            label O2GlobalId_;
            //- CO positions in global list
            label COGlobalId_;
            //- CO2 position in global list
            label CO2GlobalId_;
            //- H2O positions in global list
            label H2OGlobalId_;
            //- H2 positions in global list
            label H2GlobalId_;



        // Local copies of thermo properties

            //- Molecular weight of C [kg/kmol]
            scalar WC_;
            //- Molecular weight of O2 [kg/kmol]
            scalar WO2_;
            //- Molecular weight of CO [kg/kmol]
            scalar WCO_;
            //- Molecular weight of CO2 [kg/kmol]
            scalar WCO2_;
            //- Molecular weight of H2O [kg/kmol]
            scalar WH2O_;
            //- Molecular weight of H2 [kg/kmol]
            scalar WH2_;

            //- Formation enthalpy for CO [J/kg]
            scalar HcCO_;
            //- Formation enthalpy for CO2 [J/kg]
            scalar HcCO2_;
            //- Formation enthalpy for H2O [J/kg]
            scalar HcH2O_;

public:

    //- Runtime type information
    TypeName("CMultipleSurfaceReactions");


    // Constructors

        //- Construct from dictionary
        CMultipleSurfaceReactions
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Construct copy
        CMultipleSurfaceReactions
        (
            const CMultipleSurfaceReactions<CloudType>& srm
        );

        //- Construct and return a clone
        virtual autoPtr<SurfaceReactionModel<CloudType>> clone() const
        {
            return autoPtr<SurfaceReactionModel<CloudType>>
            (
                new CMultipleSurfaceReactions<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~CMultipleSurfaceReactions();


    // Member Functions

        //- Update surface reactions
        virtual scalar calculate
        (
            const scalar dt,
            const label celli,
            const scalar d,
            const scalar T,
            const scalar Tc,
            const scalar pc,
            const scalar rhoc,
            const scalar mass,
            const scalarField& YGas,
            const scalarField& YLiquid,
            const scalarField& YSolid,
            const scalarField& YMixture,
            const scalar N,
            scalarField& dMassGas,
            scalarField& dMassLiquid,
            scalarField& dMassSolid,
            scalarField& dMassSRCarrier
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CMultipleSurfaceReactions.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
