#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "problems.h"

using namespace std;

const double t0 = 0.0;
const double T = 100;
const double eps = 1e-6;
const double tauStart = 1;
const vector<double> u0 = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
const vector<double> lambdas = { -1.0, -0.1, -2.0 }; // это для Далквиста
string funcType = "heat";

double norm2OfVector(const vector<double>& vector){
    if (!vector.size())
        return NAN;
    double sum = 0;
    for (int i = 0; i < vector.size(); i++){
        sum += pow(vector[i], 2);
    }
    return sqrt(sum);
}

double normInfOfVector(const vector<double>& vector) {
    if (!vector.size())
        return NAN;
    double max = abs(vector[0]);
    for (int i = 1; i < vector.size(); i++)
        if (abs(vector[i]) > max)
            max = abs(vector[i]);
    return max;
}

vector<double> f(const vector<double>& u){
    if (funcType == "dahlquist") {
        return dahlquistRightSide(u, lambdas);
    }
    if (funcType == "heat") {
        return heatRightSide(u);
    }
    return {0, 0};
}

// Одна итерация алгоритма А (Новиков с.94)
void iterationA(vector<double>(*f)(const vector<double>& u), vector<double>& k1,
    vector<double>& k2, vector<double>& k3, const vector<double>& tempY, double tau, double eps,
    double g, double d, double q, double& tauNew1, double& tauNew2, vector<double>& nextY) {
    temp:
    k1 = tau * f(tempY);
    k2 = tau * f(tempY + 2.0 / 3 * k1);
    double A1 = abs(1 - 6 * g) / 4 * norm2OfVector(k2 - k1);
    double s = 0.5 * logbase(eps / A1, q);
    if (s < 0) {
        tau = pow(q, s) * tau;
        goto temp;
    }
    k3 = tau * f(tempY + 1.0 / 3 * (k1 + k2));
    nextY = tempY + (1.0 / 4) * k1 + (3 - 18.0 * g) / 4 * k2 + 9.0 * g / 2 * k3;
    double A2 = abs(1 - 6 * g) / 6 * norm2OfVector(tau * f(nextY) - k1);
    double v = 0.5 * logbase(eps / A2, q);
    if (v < 0) {
        tau = pow(q, v) * tau;
        goto temp;
    }
    double V = 3 * normInfOfVector((k3 - k2) / (k2 - k1));
    double s1 = 0.5 * logbase(eps / (d * A1), q);
    double v1 = 0.5 * logbase(eps / (d * A2), q);
    double r1 = logbase(18.0 / V, q);
    double r = logbase(6.0 / V, q);
    tauNew1 = max(tau, pow(q, min(s1, v1, r1)) * tau);
    tauNew2 = max(tau, pow(q, min(s, v, r)) * tau);
}

// Одна итерация алгоритма Б (Новиков с.95)
void iterationB(vector<double>(*f)(const vector<double>& U), vector<double>& k1,
        vector<double>& k2, vector<double>& k3, const vector<double>& tempY, double tau, double eps,
        double g, double d, double q, double& tauNew1, double& tauNew2, vector<double>& nextY) {
    temp: 
    k1 = tau * f(tempY);
    k2 = tau * f(tempY + 2.0 / 3 * k1);
    double A1 = abs(1 - 6 * g) / 4 * norm2OfVector(k2 - k1);
    double s = 0.5 * logbase(eps / (d * A1), q);
    if (s < 0) {
        tau = pow(q, s) * tau;
        goto temp;
    }
    k3 = tau * f(tempY + 1.0 / 3 * (k1 + k2));
    nextY = tempY + (7.0 / 9) * k1 + (16.0 / 81) * k2 + (2.0 / 81) * k3;
    double A2 = abs(1 - 6 * g) / 6 * norm2OfVector(tau * f(nextY) - k1);
    double v = 0.5 * logbase(eps / (d * A2), q);
    if (v < 0) {
        tau = pow(q, v) * tau;
        goto temp;
    }
    double V = 3 * normInfOfVector((k3 - k2) / (k2 - k1));
    double r1 = logbase(18.0 / V, q);
    double r = logbase(6.0 / V, q);
    tauNew1 = max(tau, pow(q, min(s, v, r1)) * tau);
    tauNew2 = max(tau, pow(q, min(s, v, r)) * tau);
}

// Метод Рунге - Кутты третьего порядка (Новиков страница с 89, 
// для него область устойчивости Q = 1 + z + 4/27 z^2 + 4/729 z^3)
void rungeKuttaMethod3WithKoef(vector<double>(*f)(const vector<double>& u), double t0, 
        double T, double tauStart, const vector<double>& u0, vector<vector<double>>& solution, 
        vector<double> &timeStamps, double eps) {
	size_t n = u0.size();
	solution.resize(1); // Здесь будут векторы-решения на каждом шаге (то есть это матрица решений)
    solution[0] = u0; 
    timeStamps.resize(1); // Здесь будут храниться моменты времени (узлы временнОй сетки)
    timeStamps[0] = t0;
    double tauNew1, tauNew2;

	vector<double> k1(n);
	vector<double> k2(n);
	vector<double> k3(n);
	double tau = tauStart; // Текущий шаг
	vector<double> tempY = u0; // Текущий y
	vector<double> nextY; // Следующий y
    const double g = 1.0 / 16, d = 152.0 / 45, q = exp(1); // Параметры метода
	double tempT = t0;
	double nextT;
    bool A = true; // Итерация по алгоритму Б или по А. Сначала начинаем с А
	while (tempT < T) {
        if (A) iterationA(f, k1, k2, k3, tempY, tau, eps, g, d, q, tauNew1, tauNew2, nextY);
        else   iterationB(f, k1, k2, k3, tempY, tau, eps, g, d, q, tauNew1, tauNew2, nextY);
        if (tauNew1 > tauNew2) {
            tau = tauNew1;
            //A = false; // Следующая итерация будет по алгоритму Б 
        } else {
            tau = tauNew2;
            A = true; // Следующая итерация будет по алгоритму А 
        }
        nextT = tempT + tau;
		timeStamps.push_back(nextT);
        solution.push_back(nextY);
		tempY = nextY;
		tempT = nextT;
	}
}

void printSolution(const vector<vector<double>>& solution) {
    ofstream myfile;
    myfile.open(funcType + "Solution.txt");
    size_t n = solution.size();
    for (size_t i = 0; i < n; ++i) {
        myfile << solution[i] << endl;
    }
    myfile.close();
}

void printTimeStamps(const vector<double>& timeStamps) {
    ofstream myfile;
    myfile.open("timeStamps.txt");
    size_t n = timeStamps.size();
    for (size_t i = 0; i < n; ++i) {
        myfile << timeStamps[i] << endl;
    }
    myfile.close();
}

void printExactSolution(const vector<double>& timeStamps) { // это только для Далквиста
    ofstream myfile;
    myfile.open(funcType + "ExactSolution.txt");
    if (funcType == "dahlquist") {
        for (size_t i = 0; i < timeStamps.size(); ++i) {
            double t = timeStamps[i];
            myfile << dahlquistExactSolution(t, u0, lambdas) << endl;
        }
    }
    myfile.close();
}

int main() {
    vector<vector<double>> solution; // Матрица решений
    vector<double> timeStamps; // Вектор из моментов времени (узлов временнОй сетки)
    rungeKuttaMethod3WithKoef(f, t0, T, tauStart, u0, solution, timeStamps, eps);
    printTimeStamps(timeStamps);
    printSolution(solution);
    
    /*printExactSolution(timeStamps);
    size_t numberOfSteps = timeStamps.size();
    size_t n = u0.size();
    vector<double> errors(numberOfSteps);
	for (size_t i = 0; i < numberOfSteps; ++i) {
		double t = timeStamps[i];
        vector<double> exactSolution = dahlquistExactSolution(t, u0, lambdas);
        errors[i] = norm2OfVector(exactSolution - solution[i]);
	}
    cout << norm2OfVector(errors) << endl;*/
    
}
