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
    casadi::MX pv_power = pvPowerModel(irradiation(time)[0]);
    std::cout << pv_power << std::endl;

	return 0;
}