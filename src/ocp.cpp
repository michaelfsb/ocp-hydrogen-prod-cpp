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
    

/*
    OPTIMAL CONTROL
*/
int main(){

    std::cout << "----->>>>> INICIO !!!!!!!!!!!!!!";


    std::vector<double> irradiationData = getDataFromFile("irradiation.txt");
    std::vector<std::vector<double>> timeData = createTimeData(irradiationData.size());
    casadi::Function irradiation = casadi::interpolant("irradiation", "bspline", timeData, irradiationData);

    std::vector<double> hydrogemDemandData = getDataFromFile("hydrogen_demand.txt");
    timeData = createTimeData(hydrogemDemandData.size());
    casadi::Function hydrogemDemand = casadi::interpolant("hydrogemDemand", "bspline", timeData, hydrogemDemandData);


    double test_pv_model = pv_model(1);

	return 0;
}