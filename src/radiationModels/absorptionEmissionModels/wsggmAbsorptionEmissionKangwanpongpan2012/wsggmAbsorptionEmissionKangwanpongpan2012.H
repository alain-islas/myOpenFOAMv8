/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::radiation::wsggmAbsorptionEmissionKangwanpongpan2012

Description
    Weighted-sum-of-gray-gases-model (WSGGM) for gas phase absorption/emission
    valid for variable (H2O/CO2) molar ratios from 0.125 to 4.

Reference
    Kangwanpongpan, T., França, F. H., da Silva, R. C., Schneider, P. S., &
    Krautz, H. J. (2012). New correlations for the weighted-sum-of-gray-gases
    model in oxy-fuel conditions based on HITEMP 2010 database. International
    Journal of Heat and Mass Transfer, 55(25-26), 7419-7433.

SourceFiles
    wsggmAbsorptionEmissionKangwanpongpan2012.C

\*---------------------------------------------------------------------------*/

#ifndef radiationwsggmAbsorptionEmissionKangwanpongpan2012_H
#define radiationwsggmAbsorptionEmissionKangwanpongpan2012_H

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
                  Class wsggmAbsorptionEmissionKangwanpongpan2012 Declaration
\*---------------------------------------------------------------------------*/

class wsggmAbsorptionEmissionKangwanpongpan2012
:
    public absorptionEmissionModel
{

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- SLG thermo package
        const fluidThermo& thermo_;

        //- Matrix coefficiens for the absorption coefficients
        scalarListList CK_;

        //- Matrix coefficients for the fiting factors
        scalarListList C_;

        //characteristic path length
        dimensionedScalar pathLength_;

        //Flag for meanBeamPath calculation
        bool meanBeamPathAutoCalcMode_;

        //Domain sector //[degree]
        dimensionedScalar sector_;



public:

    //- Runtime type information
    TypeName("wsggmAbsorptionEmissionKangwanpongpan2012");


    // Constructors

        //- Construct from components
        wsggmAbsorptionEmissionKangwanpongpan2012
	(
		const dictionary& dict,
		const fvMesh& mesh,
		const word& modelName=typeName
	);


    //- Destructor
    virtual ~wsggmAbsorptionEmissionKangwanpongpan2012();


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
