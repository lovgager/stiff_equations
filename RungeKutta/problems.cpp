#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "problems.h"


using namespace std;

// Тут переопределяем операторы

ostream& operator<<(ostream &os, const vector<vector<double>> &matrix){
    int rows = matrix.size();
    int cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else{
        os << 0;
        return os;
    }
    for (int i = 0; i < rows - 1; i++){
        for (int j = 0; j < cols; j++){
            os << matrix[i][j] << ' ';
        }
        os << '\n';
    }
    for (int j = 0; j < cols; j++){
            os << matrix[rows - 1][j] << ' ';
        }
    return os;
}

ostream& operator<<(ostream &os, const vector<double> &vector){
    int rows = vector.size();
    if (!rows){
        os << 0;
        return os;
    }  
    // os << "{ ";
    for (int i = 0; i < rows - 1; i++)
        os << vector[i] << ", ";
    os << vector[rows - 1] << ' ';
    // os << '}';
    return os;
}

vector<double> operator+(const vector<double>& vec1, const vector<double>& vec2)
{
    auto result = vec1;
    const int size_max = max<int>(vec1.size(), vec2.size());
    result.resize(size_max, 0);
    for (int i = 0; i < result.size() && i < vec2.size(); i++){   
        result[i] += vec2[i];
    }
    return result;
}

vector<double> operator-(const vector<double>& vec1, const vector<double>& vec2)
{
    auto result = vec1;
    const int size_max = max<int>(vec1.size(), vec2.size());
    result.resize(size_max, 0);
    for (int i = 0; i < result.size() && i < vec2.size(); i++){   
        result[i] -= vec2[i];
    }
    return result;
}

vector<double> operator/(const vector<double>& vec1, const vector<double>& vec2)
{
    auto result = vec1;
    const int size_max = max<int>(vec1.size(), vec2.size());
    result.resize(size_max, 0);
    for (int i = 0; i < result.size() && i < vec2.size(); i++) {
        result[i] /= vec2[i];
    }
    return result;
}

vector<double> operator*(const vector<vector<double>> &matrix, const vector<double> &vec)
{
    int rows1 = matrix.size();
    int cols = 0;
    if (rows1 != 0)
        cols = matrix[0].size();
    else
        return vector<double>(1, NAN);
    int rows2 = vec.size();
    if (cols != rows2)
        return std::vector<double>(rows2, NAN);
    vector<double> result(rows1);
    for (int i = 0; i < rows1; i++){
        double sum = 0;
        for (int k = 0; k < cols; k++){
            sum += matrix[i][k] * vec[k];
        }
        result[i] = sum;
    }
    return result;
}

vector<double> operator*(double num, const vector<double> &vec){
    int size = vec.size();
    if (!size)
        return vector<double>(1, NAN);
    vector<double> res(size);
    for (int i = 0; i < size; i++){
        res[i] = num * vec[i];
    }
    return res;
}

vector<vector<double>> operator*(double num, const vector<vector<double>> &matrix){
    int rows = matrix.size();
    int cols = 0;
    if (!rows)
        return vector<vector<double>>(1, vector<double>(1, NAN));
    else
        cols = matrix[0].size();
    vector<vector<double>> res;
    vector<double> tempVec(cols);
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++)
        {
            tempVec[j] = num * matrix[i][j];
        }
        res.push_back(tempVec);
    }
    return res;
}

double max(double a, double b) {
    return a > b ? a : b;
}

double min(double a, double b, double c) {
    double res = a;
    if (b < res) res = b;
    if (c < res) res = c;
    return res;
}

double logbase(double a, double base) {
    return log(a) / log(base);
}

vector<double> dahlquistRightSide(const vector<double> &u, const vector<double> &lambdas) {
    size_t n = u.size();
    vector<double> rightSide(n);
    for (size_t i = 0; i < n; ++i)
        rightSide[i] = lambdas[i] * u[i];
    return rightSide;
}

vector<double> dahlquistExactSolution(double t, const vector<double>& u0, const vector<double>& lambdas) {
    size_t n = lambdas.size();
    vector<double> result(n);
    for (size_t i = 0; i < n; ++i)
        result[i] = u0[i] * exp(lambdas[i] * t);
    return result;
}

vector<double> heatRightSide(const vector<double>& u) {
    size_t n = u.size();
    vector<double> rightSide(n);
    rightSide[0] = u[1] - 2 * u[0];
    for (size_t i = 1; i < n - 1; ++i)
        rightSide[i] = u[i + 1] - 2 * u[i] + u[i - 1];
    rightSide[n - 1] = -2 * u[n - 1] + u[n - 2];
    return rightSide;
}
