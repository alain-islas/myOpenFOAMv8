/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

Class
    Foam::radiation::wsggmAbsorptionEmissionSmith1982

Description
    Weighted-sum-of-gray-gases-model (WSGGM) for gas phase absorption/emission
    based on model coefficients of Smith et al., 1982
    
    Reference:                                                                          
    \verbatim                                                                           
        Smith, T. F., Z. F. Shen, and J. N. Friedman. "Evaluation of 
        coefficients for the weighted sum of gray gases model." (1982):
        602-608.                                          
    \endverbatim 
    
SourceFiles
    wsggmAbsorptionEmissionSmith1982.C

\*---------------------------------------------------------------------------*/

#ifndef radiationwsggmAbsorptionEmissionSmith1982_H
#define radiationwsggmAbsorptionEmissionSmith1982_H

#include "absorptionEmissionModel.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
/*---------------------------------------------------------------------------*\
                  Class wsggmAbsorptionEmissionSmith1982 Declaration
\*---------------------------------------------------------------------------*/

class wsggmAbsorptionEmissionSmith1982
:
    public absorptionEmissionModel
{

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- SLG thermo package
        const fluidThermo& thermo_;

        //- absorption coefficitent for individual grey gases
        scalarList emissivityCoeffs_;

        //- fitting Factors
        scalarListList fittingFactors_;

        //characteristic path length
        dimensionedScalar pathLength_;

        //Flag for meanBeamPath calculation
        bool meanBeamPathAutoCalcMode_;

        //Domain sector //[degree]
        dimensionedScalar sector_;



public:

    //- Runtime type information
    TypeName("wsggmAbsorptionEmissionSmith1982");


    // Constructors

        //- Construct from components
        wsggmAbsorptionEmissionSmith1982
	(
		const dictionary& dict,
		const fvMesh& mesh,
		const word& modelName=typeName
	);


    //- Destructor
    virtual ~wsggmAbsorptionEmissionSmith1982();


    // Member Functions

        // Access

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> aCont(const label bandI = 0) const;

                //- Emission coefficient for continuous phase
                tmp<volScalarField> eCont(const label bandI = 0) const;

                //- Emission contribution for continuous phase
                tmp<volScalarField> ECont(const label bandI = 0) const;

                //- wF for continuous phase
                tmp<volScalarField> wFCont(const label bandI = 0) const;

    // Member Functions

        inline bool isGrey() const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorptionEmissionModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
