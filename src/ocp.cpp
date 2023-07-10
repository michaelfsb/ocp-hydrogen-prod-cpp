#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <casadi/casadi.hpp>

using namespace std;
using namespace casadi;

/*
    INPUT DATA
*/
vector<vector<double>> createTimeData(const int size) {
    vector<vector<double>> timeGrid;
    vector<double> times;

    for (double time = 0.0; time < size; time++) {
        times.push_back({ time });
    }

    timeGrid.push_back(times);

    return timeGrid;
}

vector<double> getDataFromFile(const string& filename) {
    vector<double> dataSerie;
    ifstream file(filename);
    if (!file) {
        cerr << "Failed to open file: " << filename << endl;
        return dataSerie;
    }

    string line;

    while (getline(file, line)) {
        try {
            double value = stod(line);
            dataSerie.push_back(value);
        }
        catch (const exception&) {
            // Ignore the line if the conversion to double fails
        }
    }

    file.close();

    return dataSerie;
}

vector<casadi_int> createIntVector(double start, double end, double step) {
    vector<casadi_int> numbers;
    for (double value = start; value <= end; value += step) {
        numbers.push_back(value);
    }
    return numbers;
}

/*
    PHOTOVOLTAIC PANEL MODEL
*/
MX lambertw(const MX& x) {
    const double E = 0.4586887;
	return (1 + E) * log(6. / 5 * x / log(12. / 5 * x / log(1 + 12. / 5 * x))) - E * log(2 * x / log(1 + 2 * x));
}

MX pvPowerModel(const MX& irradiation) {
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
    MX Iph = (Isc + Kl * (T_ps - Tr)) * irradiation;

    // Algebraic equations
    MX v_ps = (N_ss * Vt * A * (lambertw(exp(1) * (Iph / Irs + 1)) - 1));
    MX i_ps = N_ps * (Iph - Irs * (exp(v_ps / (N_ss * Vt * A)) - 1));

   return i_ps * v_ps;
}

MX electrolyzerProdModel(const MX& i_el) {
    // Declare constants
    double F = 96485.33289;  // Faraday constant

    // Declare electrolyzer parameters
    int N_el = 6;            // Number of cells

    return (N_el * i_el / (F * 1000)) * 60 * 11.126; // Hydrogen production rate
}

/*
    PHOTOVOLTAIC PANEL MODEL
*/
MX electrolyzerPowerModel(const MX& i_el) {
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
    MX i_el_d = i_el / A_el; // Current density
    double ro_b = (0.005139 * lambda_b - 0.00326) * exp(1268 * (1 / 303 - 1 / t_el)); // Membrane conductivity
    double v_el_0 = 1.23 - 0.0009 * (t_el - 298) + 2.3 * R * t_el * log(P_h2 * P_h2 * P_o2) / (4 * F); // Reversible potential of the electrolyzer
    MX v_etd = (R * t_el / F) * asinh(0.5 * i_el_d / I_ao) + (R * t_el / F) * asinh(0.5 * i_el_d / I_co) + i_el_d * delta_b / ro_b; // Electrode overpotential
    MX v_el_hom_ion = delta_b * i_el / (A_el * ro_b); // Ohmic overvoltage and ionic overpotential

    MX v_el = N_el * (v_el_0 + v_etd + v_el_hom_ion); // Electrolyzer voltage
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
    MX v_h2 = MX::sym("v_h2"); // State - Volume of hydrogen
    MX i_el = MX::sym("i_el"); // Control - Electrical current in electrolyzer
    MX time = MX::sym("time"); // Time

    // Input data
    vector<double> irradiationData = getDataFromFile("irradiation.txt");
    vector<vector<double>> timeData = createTimeData(irradiationData.size());
    Function irradiation = interpolant("irradiation", "bspline", timeData, irradiationData);

    vector<double> hydrogemDemandData = getDataFromFile("hydrogen_demand.txt");
    timeData = createTimeData(hydrogemDemandData.size());
    Function hydrogemDemand = interpolant("hydrogemDemand", "bspline", timeData, hydrogemDemandData);

    // Models
    MX p_pv = pvPowerModel(irradiation(time)[0]);
    MX p_el = electrolyzerPowerModel(i_el);
    MX f_h2 = electrolyzerProdModel(i_el);
    MX v_h2_dot = 1 - hydrogemDemand(time)[0];

    // Lagrange cost function
    MX f_l = pow((p_el - p_pv), 2);

    // Creat NPL problem
    double h = Tf / N;
    
    vector<double> t(N);
    for (int i = 0; i < N; ++i) {
        t[i] = i * (Tf / (N - 1));
    }

    Function f = Function("f", 
        { v_h2, i_el, time }, { v_h2_dot, f_l }, 
        { "x", "u", "t" }, { "x_dot", "L" });

    MXVector X, U;
    for (double k = 0; k < N - 0.5; k += 0.5) {
        X.push_back(MX::sym("X_" + to_string(k)));
        U.push_back(MX::sym("U_" + to_string(k)));
    }

    MX L = MX(0.0);
    MXVector g;
    MXVector lbg;
    MXVector ubg;
    MXVector w;
    MXVector w0;
    MXVector lbw;
    MXVector ubw;

    for (double k = 0; k < 2 * N - 2; k += 2) {
        int i = static_cast<int>(k / 2);

        MXVector fw0 = f({ X[k], U[k], t[i] });
        MXVector fw1 = f({ X[k + 1], U[k + 1], t[i] + h / 2 });
        MXVector fw2 = f({ X[k + 2], U[k + 2], t[i + 1] });

        L = L + h * (fw0[1] + 4 * fw1[1] + fw2[1]) / 6;

        // Add equality constraint
        g.push_back(X[k + 2] - X[k] - h * (fw0[0] + 4 * fw1[0] + fw2[0]) / 6);
        lbg.push_back(MX(0));
        ubg.push_back(MX(0));
        g.push_back(X[k + 1] - (X[k + 2] + X[k]) / 2 - h * (fw0[0] - fw2[0]) / 8);
        lbg.push_back(MX(0));
        ubg.push_back(MX(0));
    }

    for (int k = 0; k < 2 * N - 1; ++k) {
        w.push_back(U[k]);
        lbw.push_back(MX(I_e_min));
        ubw.push_back(MX(I_e_max));
        w0.push_back(MX(I_e_0));

        w.push_back(X[k]);
        lbw.push_back(MX(M_min));
        ubw.push_back(MX(M_max));
        w0.push_back(MX(M_0));
    }

    // Set the initial condition for the state
    lbw[1] = M_0;
    ubw[1] = M_0;

    MX w_concat = w[0];
    for (int i = 1; i < w.size(); ++i) {
        w_concat = horzcat(w_concat, w[i]);
    }

    MX g_concat = g[0];
    for (int i = 1; i < g.size(); ++i) {
        g_concat = horzcat(g_concat, g[i]);
    }

    Function solver = nlpsol("solver", "ipopt", { {"f", L}, {"x", w_concat}, {"g", g_concat} });

    // Solve the problem
    DMDict arg = { {"lbx", lbw},
                       {"ubx", ubw},
                       {"x0", w0},
                       {"lbg", lbg},
                       {"ubg", ubg} };
    DMDict res = solver(arg);

    // Print the optimal cost
    double cost(res.at("f"));
    std::cout << "optimal cost: " << cost << std::endl;

    // Print the optimal solution
    std::vector<double> uopt(res.at("x"));
    std::cout << "optimal control: " << uopt << std::endl;

	return 0;
}