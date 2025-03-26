#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>

using namespace std;

void Check(const vector<double>& alpha, const vector<double>& beta, const vector<double>& y) {
    string res = "Достаточные условия выполнены!!!";
    for (size_t i = 0; i < y.size() - 1; ++i) {
        if (abs(y[i]) < abs(alpha[i]) + abs(beta[i])) {
            res = "Достаточные условия не выполнены!!!";
        }
    }
    cout << res << endl;
}

void printMatrixA(const vector<double>& a, const vector<double>& b, const vector<double>& c) {
    int n = b.size();
    cout << "Матрица A:" << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                cout << setw(4) << b[i] << " ";
            } else if (j == i + 1 && i < n - 1) {
                cout << setw(4) << c[i] << " ";
            } else if (i == j + 1 && j < n - 1) {
                cout << setw(4) << a[i] << " ";
            } else {
                cout << setw(4) << 0 << " ";
            }
        }
        cout << endl;
    }
}

string scientificFormat(double value, int precision = 14) {
    ostringstream oss;
    oss << scientific << setprecision(precision) << value;
    string s = oss.str();

    // Удаляем лишние нули в экспоненте
    size_t e_pos = s.find('e');
    if (e_pos != string::npos) {
        size_t last_non_zero = s.find_last_not_of('0', e_pos - 1);
        if (last_non_zero != string::npos && s[last_non_zero] != '.') {
            s.erase(last_non_zero + 1, e_pos - last_non_zero - 1);
        }
    }
    return s;
}

void printVectorColumn(const string& title, const vector<double>& vec, bool scientific_fmt = false) {
    cout << title << endl;
    for (size_t i = 0; i < vec.size(); ++i) {
        if (scientific_fmt) {
            cout <<  scientificFormat(vec[i]) << endl;
        } else {
            cout << fixed << setprecision(14) << vec[i] << endl;
        }
    }
    cout << endl;
}

int main() {
    // Исходные данные
    vector<double> a = {2, 1, 0, 1, 1, 2};
    vector<double> b = {95, 54, 81, 131, 67, 91, 84};
    vector<double> c = {-1, 1, 0, 2, -3, 0};
    vector<double> d = {17, 9, 15, 24, 13, 17, 19};

    // Добавляем элементы в векторы a и c
    a.insert(a.begin(), 0);
    c.push_back(0);

    // Выводим матрицу A
    printMatrixA(a, b, c);
    cout << endl;

    // Прямая прогонка
    vector<double> alpha(7, 0);
    vector<double> beta(7, 0);
    vector<double> y(7, 0);

    y[0] = b[0];
    alpha[0] = -c[0] / y[0];
    beta[0] = d[0] / y[0];

    for (int i = 1; i < 7; ++i) {
        y[i] = b[i] + a[i] * alpha[i - 1];
        alpha[i] = -c[i] / y[i];
        beta[i] = (d[i] - a[i] * beta[i - 1]) / y[i];
    }

    Check(alpha, beta, y);
    cout << endl;

    // Обратная прогонка
    vector<double> x(7, 0);
    x[6] = beta[6];
    for (int i = 1; i < 7; ++i) {
        x[6 - i] = alpha[6 - i] * x[6 - i + 1] + beta[6 - i];
    }

    // Вывод значений
    cout << "Альфа коэффициент" << endl;
    for (double val : alpha) {
        cout << fixed << setprecision(3) << setw(7) << val << " ";
    }
    cout << endl << endl;

    cout << "Бета коэффициент" << endl;
    for (double val : beta) {
        cout << fixed << setprecision(3) << setw(7) << val << " ";
    }
    cout << endl << endl;

    // Вывод вектора x в столбик
    printVectorColumn("Решение x:", x,false);

    // Расчет невязки
    vector<double> r(7, 0);
    vector<vector<double>> A = {
        {95, -1, 0, 0, 0, 0, 0},
        {2, 54, 1, 0, 0, 0, 0},
        {0, 1, 81, 0, 0, 0, 0},
        {0, 0, 0, 131, 2, 0, 0},
        {0, 0, 0, 1, 67, -3, 0},
        {0, 0, 0, 0, 1, 91, 0},
        {0, 0, 0, 0, 0, 2, 84}
    };

    for (int i = 0; i < 7; ++i) {
        double sum = 0;
        for (int j = 0; j < 7; ++j) {
            sum += A[i][j] * x[j];
        }
        r[i] = d[i] - sum;
    }

    // Вывод невязки в научном формате
    printVectorColumn("Невязка СЛАУ:", r, true);

    // Нормы невязки
    double norm_1 = 0;
    double norm_inf = 0;
    for (double val : r) {
        norm_1 += abs(val);
        if (abs(val) > norm_inf) {
            norm_inf = abs(val);
        }
    }

    cout << fixed << setprecision(14);
    cout << "Единичная норма невязки: " << norm_1 << endl;
    cout << "Бесконечная норма невязки: " << norm_inf<< endl << endl;

    // Устойчивость к малым возмущениям
    vector<double> d1 = d;
    for (double& val : d1) {
        val += 0.01;
    }

    // Повторяем прогонку с измененным вектором d1
    vector<double> alpha1(7, 0);
    vector<double> beta1(7, 0);
    vector<double> y1(7, 0);

    y1[0] = b[0];
    alpha1[0] = -c[0] / y1[0];
    beta1[0] = d1[0] / y1[0];

    for (int i = 1; i < 7; ++i) {
        y1[i] = b[i] + a[i] * alpha1[i - 1];
        alpha1[i] = -c[i] / y1[i];
        beta1[i] = (d1[i] - a[i] * beta1[i - 1]) / y1[i];
    }

    vector<double> x1(7, 0);
    x1[6] = beta1[6];
    for (int i = 1; i < 7; ++i) {
        x1[6 - i] = alpha1[6 - i] * x1[6 - i + 1] + beta1[6 - i];
    }

    cout << "Устойчивость к малым возмущениям" << endl;
    cout << "Альфа коэффициент" << endl;
    for (double val : alpha1) {
        cout << fixed << setprecision(3) << setw(7) << val << " ";
    }
    cout << endl << endl;

    cout << "Бета коэффициент" << endl;
    for (double val : beta1) {
        cout << fixed << setprecision(3) << setw(7) << val << " ";
    }
    cout << endl << endl;

    // Вывод возмущенного решения в столбик
    printVectorColumn("Возмущенное решение x1:", x1, false);

    return 0;
}