#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <casadi/casadi.hpp>

/*
    INPUT DATA
*/
std::vector<std::vector<double>> createTimeData(const int size) {
    std::vector<std::vector<double>> timeGrid;
    std::vector<double> times;

    for (double time = 0.0; time < size; time++) {
        times.push_back({ time });
    }

    timeGrid.push_back(times);

    return timeGrid;
}

std::vector<double> getDataFromFile(const std::string& filename) {
    std::vector<double> dataSerie;
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return dataSerie;
    }

    std::string line;

    while (std::getline(file, line)) {
        try {
            double value = std::stod(line);
            dataSerie.push_back(value);
        }
        catch (const std::exception&) {
            // Ignore the line if the conversion to double fails
        }
    }

    file.close();

    return dataSerie;
}

std::vector<casadi_int> createIntVector(double start, double end, double step) {
    std::vector<casadi_int> numbers;
    for (double value = start; value <= end; value += step) {
        numbers.push_back(value);
    }
    return numbers;
}

/*
    PHOTOVOLTAIC PANEL MODEL
*/
casadi::MX lambertw(const casadi::MX& x) {
    const double E = 0.4586887;
	return (1 + E) * log(6. / 5 * x / log(12. / 5 * x / log(1 + 12. / 5 * x))) - E * log(2 * x / log(1 + 2 * x));
}

casadi::MX pvPowerModel(const casadi::MX& irradiation) {
    // Constants
    const double Q = 1.6e-19; // Elementary charge
    const double K = 1.38e-23; // Boltzmann constant

    // Photovoltaic parameters
    const double N_ps = 8;           // Number of panels in parallel
    const double N_ss = 300;         // Number of panels in series
    const double T_ps = 298;         // Temperature
    const double Tr = 298;           // Reference temperature
    const double Isc = 3.27;         // Short circuit current at Tr
    const double Kl = 0.0017;        // Short circuit current temperature coeff
    const double Ior = 2.0793e-6;    // Ior - Irs at Tr
    const double Ego = 1.1;          // Band gap energy of the semiconductor
    const double A = 1.6;            // Factor.cell deviation from de ideal pn junction

    // Intermediate variables
    const double Vt = K * T_ps / Q;
    const double Irs = Ior * pow((T_ps / Tr),3) * exp(Q * Ego * (1 / Tr - 1 / T_ps) / (K * A));
    casadi::MX Iph = (Isc + Kl * (T_ps - Tr)) * irradiation;

    // Algebraic equations
    casadi::MX v_ps = (N_ss * Vt * A * (lambertw(exp(1) * (Iph / Irs + 1)) - 1));
    casadi::MX i_ps = N_ps * (Iph - Irs * (exp(v_ps / (N_ss * Vt * A)) - 1));

   return i_ps * v_ps;
}

casadi::MX electrolyzerProdModel(const casadi::MX& i_el) {
    // Declare constants
    double F = 96485.33289;  // Faraday constant

    // Declare electrolyzer parameters
    int N_el = 6;            // Number of cells

    return (N_el * i_el / (F * 1000)) * 60 * 11.126; // Hydrogen production rate
}

/*
    PHOTOVOLTAIC PANEL MODEL
*/
casadi::MX electrolyzerPowerModel(const casadi::MX& i_el) {
    // Declare constants
    const double R = 8.314;                 // Gas constant
    const double F = 96485.33289;           // Faraday constant

    // Declare electrolyzer parameters
    const double A_el = 212.5;              // Stack area
    const int N_el = 6;                     // Number of cells
    const double P_h2 = 6.9;                // Hydrogen partial pressure
    const double P_o2 = 1.3;                // Oxygen partial pressure
    const double I_ao = 1.0631e-6;          // Anode current density
    const double I_co = 1e-3;               // Cathode current density
    const double delta_b = 178e-6;          // Membrane thickness
    const int lambda_b = 21;                // Membrane water content
    const int t_el = 298;                   // Temperature

    // Intermediate electrolyzer variables
    casadi::MX i_el_d = i_el / A_el; // Current density
    double ro_b = (0.005139 * lambda_b - 0.00326) * exp(1268 * (1 / 303 - 1 / t_el)); // Membrane conductivity
    double v_el_0 = 1.23 - 0.0009 * (t_el - 298) + 2.3 * R * t_el * log(P_h2 * P_h2 * P_o2) / (4 * F); // Reversible potential of the electrolyzer
    casadi::MX v_etd = (R * t_el / F) * asinh(0.5 * i_el_d / I_ao) + (R * t_el / F) * asinh(0.5 * i_el_d / I_co) + i_el_d * delta_b / ro_b; // Electrode overpotential
    casadi::MX v_el_hom_ion = delta_b * i_el / (A_el * ro_b); // Ohmic overvoltage and ionic overpotential

    casadi::MX v_el = N_el * (v_el_0 + v_etd + v_el_hom_ion); // Electrolyzer voltage
    return i_el * v_el; // Electrolyzer consumed power
}

/*
    OPTIMAL CONTROL
*/
int main(){
    // Parameters of OCP
    const int N = 15;             // Number of control intervals
    const double Tf = 1440.0;     // Final time (min)

    // Intial and limit values
    const double M_0 = 0.65;      // Initial volume of hydrogen (Nm3)
    const double M_min = 0.6;     // Minimum volume of hydrogen (Nm3)
    const double M_max = 2.5;     // Maximum volume of hydrogen (Nm3)
    const double I_e_0 = 30.0;    // Initial current (A)
    const double I_e_min = 1.0;   // Minimum current (A)
    const double I_e_max = 100.0; // Maximum current (A)
 
    // Variables of OCP
    casadi::MX v_h2 = casadi::MX::sym("v_h2"); // State - Volume of hydrogen
    casadi::MX i_el = casadi::MX::sym("i_el"); // Control - Electrical current in electrolyzer
    casadi::MX time = casadi::MX::sym("time"); // Time

    // Input data
    std::vector<double> irradiationData = getDataFromFile("irradiation.txt");
    std::vector<std::vector<double>> timeData = createTimeData(irradiationData.size());
    casadi::Function irradiation = casadi::interpolant("irradiation", "bspline", timeData, irradiationData);

    std::vector<double> hydrogemDemandData = getDataFromFile("hydrogen_demand.txt");
    timeData = createTimeData(hydrogemDemandData.size());
    casadi::Function hydrogemDemand = casadi::interpolant("hydrogemDemand", "bspline", timeData, hydrogemDemandData);

    // Models
    casadi::MX p_pv = pvPowerModel(irradiation(time)[0]);
    casadi::MX p_el = electrolyzerPowerModel(i_el);
    casadi::MX f_h2 = electrolyzerProdModel(i_el);
    casadi::MX v_h2_dot = 1 - hydrogemDemand(time)[0];

    // Lagrange cost function
    casadi::MX f_l = pow((p_el - p_pv), 2);

	return 0;
}