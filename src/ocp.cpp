#include <iostream>
#include <casadi/casadi.hpp>

using namespace casadi;

double lambertw(const double& x) {
    const double E = 0.4586887;
	return (1 + E) * log(6. / 5 * x / log(12. / 5 * x / log(1 + 12. / 5 * x))) - E * log(2 * x / log(1 + 2 * x));
}

double pv_model(const double& irradiation) {
    // Declare constants
    const double Q = 1.6e-19; // Elementary charge
    const double K = 1.38e-23; // Boltzmann constant

    // Declare photovoltaic parameters
    const double N_ps = 8;           // Number of panels in parallel
    const double N_ss = 300;         // Number of panels in series
    const double T_ps = 298;         // Temperature
    const double Tr = 298;           // Reference temperature
    const double Isc = 3.27;         // Short circuit current at Tr
    const double Kl = 0.0017;        // Short circuit current temperature coeff
    const double Ior = 2.0793e-6;    // Ior - Irs at Tr
    const double Ego = 1.1;          // Band gap energy of the semiconductor
    const double A = 1.6;            // Factor.cell deviation from de ideal pn junction

    // Intermediate photovoltaic variables
    const double Vt = K * T_ps / Q;
    const double Irs = Ior * pow((T_ps / Tr),3) * exp(Q * Ego * (1 / Tr - 1 / T_ps) / (K * A));
    double Iph = (Isc + Kl * (T_ps - Tr)) * irradiation;


    // Algebraic equations
    double v_ps = (N_ss * Vt * A * (lambertw(exp(1) * (Iph / Irs + 1)) - 1));
    double i_ps = N_ps * (Iph - Irs * (exp(v_ps / (N_ss * Vt * A)) - 1));

   return i_ps;
}
    



int main(){
    double test = pv_model(1);

	return 0;
}