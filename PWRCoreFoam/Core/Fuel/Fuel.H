/*----------------------------------------------------------------------------*\
| Class:                                                                       |
|     Foam::Fuel                                                               |
\*----------------------------------------------------------------------------*/

#ifndef Fuel_H
#define Fuel_H

#include "fvCFD.H"
#include "rhoThermo.H"
#include "rodFieldsFwd.H"
#include "subChannelMesh.H"


//****************************************************************************//

namespace Foam
{

/*----------------------------------------------------------------------------*\
|                       Class Fuel Declaration                                 |
\*----------------------------------------------------------------------------*/

class Fuel
{

public:

    // Public typedefs


private:

    // Private data

        const IOdictionary fuelDict_;

        const subChannelMesh & mesh_;

        const volVectorField & U_;

        const rhoThermo & thermo_;

        volScalarField Cp_;

        volScalarField kappa_;

        rodScalarField rodLinePowerDensity_;

        rodScalarField pelletCenterTemperature_;

        rodScalarField cladSurfaceTemperature_;

        HashTable<label, word> powerTypeTable_;

        rodScalarField Re_for_rod_;

        rodScalarField h_for_rod_;

        rodScalarField U_time_rho_;

        rodScalarField T_coolant_;


    // Private member functions

        void set_uniform_power_distribution(scalar uniform_line_power_density);

        void set_cos_power_distribution(scalar average_line_power_density);

        void set_polynominal_power_distribution(scalar average_line_power_density);

        scalar calc_value_of_polynominal(const scalar & x,
                                         const scalarList & polynominal_coeffs);

        scalar calc_integral_of_polynominal(const scalar & x_1,
                                            const scalar & x_2,
                                            const scalarList & polynominal_coeffs);

        void set_horizontal_power_distribution();

        void set_fuel_power();

        void update_fuel_power();

        void update_clad_temperature();

        void update_pellet_temperature();

        scalar get_avg_T_around_rod(label edge_i);

        scalar get_avg_U_z_around_rod(label edge_i);

        scalar get_avg_rho_around_rod(label edge_i);

        scalar get_avg_Cp_around_rod(label edge_i);

        scalar get_avg_mu_around_rod(label edge_i);

        scalar get_avg_kappa_around_rod(label edge_i);

        scalar get_avg_alpha_around_rod(label edge_i);


public:

    // Runtime type information

        TypeName("Fuel");


    // Public typedefs


    // Static member functions


    // Constructors

        // Construct from subChannelMesh
        Fuel
        (
          const subChannelMesh & mesh,
          const volVectorField & U,
          const rhoThermo & thermo
        );

    // Destructor

        virtual ~Fuel();


    // Public member functions

        void update();

        void write();

        inline const IOdictionary & get_fuel_dict();
        inline const IOdictionary & get_fuel_dict() const;

        inline rodScalarField & get_rod_linear_power_density();
        inline const rodScalarField & get_rod_linear_power_density() const;

        inline rodScalarField & get_pellet_center_temperature();
        inline const rodScalarField & get_pellet_center_temperature() const;

        inline rodScalarField & get_clad_surface_temperature();
        inline const rodScalarField & get_clad_surface_temperature() const;


    // Member operators


    // IOstream operators


};

}  // End namespace Foam

//****************************************************************************//

#include "FuelI.H"

//****************************************************************************//

#endif


