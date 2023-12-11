/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "CMultipleSurfaceReactions.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CMultipleSurfaceReactions<CloudType>::
CMultipleSurfaceReactions
(
    const dictionary& dict,
    CloudType& owner
)
:
    SurfaceReactionModel<CloudType>(dict, owner, typeName),
    C11_(this->coeffDict().template lookup<scalar>("C11")),
    C21_(this->coeffDict().template lookup<scalar>("C21")),
    E1_(this->coeffDict().template lookup<scalar>("E1")),
    C12_(this->coeffDict().template lookup<scalar>("C12")),
    C22_(this->coeffDict().template lookup<scalar>("C22")),
    E2_(this->coeffDict().template lookup<scalar>("E2")),
    C13_(this->coeffDict().template lookup<scalar>("C13")),
    C23_(this->coeffDict().template lookup<scalar>("C23")),
    E3_(this->coeffDict().template lookup<scalar>("E3")),
    CsLocalId_(-1),
    O2GlobalId_(owner.composition().carrierId("O2")),
    COGlobalId_(owner.composition().carrierId("CO")),
    CO2GlobalId_(owner.composition().carrierId("CO2")),
    H2OGlobalId_(owner.composition().carrierId("H2O")),
    H2GlobalId_(owner.composition().carrierId("H2")),
    WC_(0.0),
    WO2_(0.0),
    WCO_(0.0),
    WCO2_(0.0),
    WH2O_(0.0),
    WH2_(0.0),
    HcCO_(0.0),
    HcCO2_(0.0),
    HcH2O_(0.0)
{
    // Determine Cs ids
    label idSolid = owner.composition().idSolid();
    CsLocalId_ = owner.composition().localId(idSolid, "C");

    // Set local copies of thermo properties
    WO2_ = owner.thermo().carrier().Wi(O2GlobalId_);
    WCO_ = owner.thermo().carrier().Wi(COGlobalId_);
    WCO2_ = owner.thermo().carrier().Wi(CO2GlobalId_);
    WH2O_ = owner.thermo().carrier().Wi(H2OGlobalId_);
    WH2_ = owner.thermo().carrier().Wi(H2GlobalId_);
    WC_ = WCO2_ - WO2_;

    HcCO_ = owner.thermo().carrier().Hf(COGlobalId_);
    HcCO2_ = owner.thermo().carrier().Hf(CO2GlobalId_);
    HcH2O_ = owner.thermo().carrier().Hf(H2OGlobalId_);

    const scalar YCloc = owner.composition().Y0(idSolid)[CsLocalId_];
    const scalar YSolidTot = owner.composition().YMixture0()[idSolid];
    Info<< "    C(s): particle mass fraction = " << YCloc*YSolidTot << endl;
}


template<class CloudType>
Foam::CMultipleSurfaceReactions<CloudType>::
CMultipleSurfaceReactions
(
    const CMultipleSurfaceReactions<CloudType>& srm
)
:
    SurfaceReactionModel<CloudType>(srm),
    C11_(srm.C11_),
    C21_(srm.C21_),
    E1_(srm.E1_),
    C12_(srm.C12_),
    C22_(srm.C22_),
    E2_(srm.E2_),
    C13_(srm.C13_),
    C23_(srm.C23_),
    E3_(srm.E3_),
    CsLocalId_(srm.CsLocalId_),
    O2GlobalId_(srm.O2GlobalId_),
    COGlobalId_(srm.COGlobalId_),
    CO2GlobalId_(srm.CO2GlobalId_),
    H2OGlobalId_(srm.H2OGlobalId_),
    H2GlobalId_(srm.H2GlobalId_),
    WC_(srm.WC_),
    WO2_(srm.WO2_),
    WCO_(srm.WCO_),
    WCO2_(srm.WCO2_),
    WH2O_(srm.WH2O_),
    WH2_(srm.WH2_),
    HcCO_(srm.HcCO_),
    HcCO2_(srm.HcCO2_),
    HcH2O_(srm.HcH2O_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CMultipleSurfaceReactions<CloudType>::
~CMultipleSurfaceReactions()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::CMultipleSurfaceReactions<CloudType>::calculate
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
) const
{
    // Fraction of remaining combustible material
    const label idSolid = CloudType::parcelType::SLD;
    const scalar fComb = YMixture[idSolid]*YSolid[CsLocalId_];

    // Surface combustion active combustible fraction is consumed
    if (fComb < small)
    {
        return 0.0;
    }

    const SLGThermo& thermo = this->owner().thermo();

    // Local mass fraction of O2 in the carrier phase
    const scalar YO2 = thermo.carrier().Y(O2GlobalId_)[celli];
   // Local mass fraction of CO2 in the carrier phase
    const scalar YCO2 = thermo.carrier().Y(CO2GlobalId_)[celli];
   // Local mass fraction of H2O in the carrier phase
    const scalar YH2O = thermo.carrier().Y(H2OGlobalId_)[celli];

    // Diffusion rate coefficient R1
    const scalar D01 = C11_/d*pow(0.5*(T + Tc), 0.75);
    // Diffusion rate coefficient R2
    const scalar D02 = C12_/d*pow(0.5*(T + Tc), 0.75);
    // Diffusion rate coefficient R3
    const scalar D03 = C13_/d*pow(0.5*(T + Tc), 0.75);

    // Kinetic rate R1
    const scalar Rk1 = C21_*exp(-E1_/(RR*Tc));
    // Kinetic rate R2
    const scalar Rk2 = C22_*exp(-E2_/(RR*Tc));
    // Kinetic rate R3
    const scalar Rk3 = C23_*exp(-E3_/(RR*Tc));

    // Particle surface area
    const scalar Ap = constant::mathematical::pi*sqr(d);

    // Change in C mass [kg]
    scalar dmC1 = Ap*rhoc*RR*Tc*YO2/WO2_*D01*Rk1/(D01 + Rk1)*dt;
    scalar dmC2 = Ap*rhoc*RR*Tc*YCO2/WCO2_*D02*Rk2/(D02 + Rk2)*dt;
    scalar dmC3 = Ap*rhoc*RR*Tc*YH2O/WH2O_*D03*Rk3/(D03 + Rk3)*dt;

    scalar dmC = dmC1 + dmC2 + dmC3;

    // Limit mass transfer by availability of C
    dmC = min(mass*fComb, dmC);

    // Molar consumption
    const scalar dOmega1 = dmC1/WC_;
    const scalar dOmega2 = dmC2/WC_;
    const scalar dOmega3 = dmC3/WC_;

    // Change in O2 mass [kg]
    const scalar dmO2 = dOmega1*WO2_;
    // Change in CO2 mass [kg] in R2
    const scalar dmCO2_R2 = dOmega2*WCO2_;
    // Change in H2O mass [kg]
    const scalar dmH2O = dOmega3*WH2O_;

    // Mass of newly created CO2 [kg] in R1
    const scalar dmCO2_R1 = dOmega1*(WCO2_);
    // Mass of newly created CO [kg] in R2
    const scalar dmCO_R2 = 2*dOmega2*(WCO_);
    // Mass of newly created CO [kg] in R3
    const scalar dmCO_R3 = dOmega3*(WCO_);
    // Mass of newly created H2 [kg]
    const scalar dmH2 = dOmega3*(WH2_);


    // Update local particle C mass
    dMassSolid[CsLocalId_] += (dOmega1 + dOmega2 + dOmega3)*WC_;

    // Update carrier consumption gases O2,CO2,H2O
    dMassSRCarrier[O2GlobalId_] -= dmO2;
    dMassSRCarrier[CO2GlobalId_] -= dmCO2_R2;
    dMassSRCarrier[H2OGlobalId_] -= dmH2O;

    // Update carrier created gases CO2,CO,H2
    dMassSRCarrier[CO2GlobalId_] += dmCO2_R1;
    dMassSRCarrier[COGlobalId_] += dmCO_R2 + dmCO_R3;
    dMassSRCarrier[H2GlobalId_] += dmH2;

    const scalar HsC = thermo.solids().properties()[CsLocalId_].Hs(T);

    // carrier sensible enthalpy exchange handled via change in mass

    // Heat of reaction [J]
    scalar Hr1 = dmC1*HsC - dmCO2_R1*HcCO2_;
    scalar Hr2 = -2*dmCO_R2*HcCO_ + dmC2*HsC + dmCO2_R2*HcCO2_;
    scalar Hr3 = -dmCO_R3*HcCO_ + dmC3*HsC + dmH2O*HcH2O_;

    scalar Hr = Hr1 + Hr2 + Hr3;
    return Hr;
}


// ************************************************************************* //
