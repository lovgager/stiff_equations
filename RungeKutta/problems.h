#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;

ostream& operator<<(ostream &os, const vector<vector<double>> &matrix);

ostream& operator<<(ostream &os, const vector<double> &vector);

vector<double> operator+(const vector<double>& vec1, const vector<double>& vec2);

vector<double> operator-(const vector<double>& vec1, const vector<double>& vec2);

vector<double> operator/(const vector<double>& vec1, const vector<double>& vec2);

vector<double> operator*(const vector<vector<double>> &matrix, const vector<double> &vec);

vector<double> operator*(double num, const vector<double> &vec);

vector<vector<double>> operator*(double num, const vector<vector<double>> &matrix);

vector<vector<double>> transpose(const vector<vector<double>>& matrix);

double max(double a, double b);

double min(double a, double b, double c);

double logbase(double a, double base);

vector<double> dahlquistRightSide(const vector<double>& u, const vector<double>& lambdas);

vector<double> dahlquistExactSolution(double t, const vector<double>& u0, const vector<double>& lambdas);

vector<double> heatRightSide(const vector<double>& u);