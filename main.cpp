#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>

using std::cout;
using std::endl;
using std::setw;
using std::vector;
using std::string;

string formatNumber(double value, int precision = 5, double pos = 1e-15) {
    std::ostringstream ss;

    if (fabs(value) < pos) {
        ss << std::fixed << std::setprecision(precision) << value;
        return ss.str();
    }

    if (fabs(value) < 1e-4 || fabs(value) >= 1e5) {
        ss << std::scientific << std::setprecision(precision) << value;
        string s = ss.str();
        auto e_pos = s.find("e+00");
        if (e_pos != string::npos || (e_pos = s.find("e-00")) != string::npos) {
            return s.substr(0, s.find("e"));
        }
        return s;
    }

    ss << std::fixed << std::setprecision(precision) << value;
    return ss.str();
}

void printVector(const vector<double> &v) {
    for (const auto &x: v)
        cout << setw(15) << formatNumber(x) << endl;
}

void printMatrix(const vector<vector<double> > &M, double pos = 1e-15) {
    for (const auto &row: M) {
        for (const auto &x: row)
            cout << setw(15) << formatNumber(x, 5, pos) << " ";
        cout << endl;
    }
}

vector<double> matVecMult(const vector<vector<double> > &A, const vector<double> &x) {
    const int n = A.size(), m = A[0].size();
    vector<double> y(n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            y[i] += A[i][j] * x[j];
    return y;
}

vector<double> vectorDiff(const vector<double> &a, const vector<double> &b) {
    const int n = a.size();
    vector<double> diff(n);
    for (int i = 0; i < n; ++i)
        diff[i] = a[i] - b[i];
    return diff;
}

double norm1_vector(const vector<double> &v) {
    double sum = 0.0;
    for (const auto &x: v)
        sum += fabs(x);
    return sum;
}

double normInf_vector(const vector<double> &v) {
    double maxVal = 0.0;
    for (const auto &x: v)
        maxVal = std::max(maxVal, fabs(x));
    return maxVal;
}

double norm1_matrix(const vector<vector<double> > &A) {
    const int n = A.size(), m = A[0].size();
    double maxSum = 0.0;
    for (int j = 0; j < m; ++j) {
        double sum = 0.0;
        for (int i = 0; i < n; ++i)
            sum += fabs(A[i][j]);
        maxSum = std::max(maxSum, sum);
    }
    return maxSum;
}

double normInf_matrix(const vector<vector<double> > &A) {
    double maxSum = 0.0;
    for (const auto &row: A) {
        double sum = 0.0;
        for (const auto &val: row)
            sum += fabs(val);
        maxSum = std::max(maxSum, sum);
    }
    return maxSum;
}

vector<vector<double> > matMult(const vector<vector<double> > &A, const vector<vector<double> > &B) {
    const int n = A.size(), m = B[0].size(), inner = A[0].size();
    vector<vector<double> > C(n, vector<double>(m, 0.0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            for (int k = 0; k < inner; ++k)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

vector<double> Gauss(const vector<vector<double> > &A, const vector<double> &F) {
    const int n = A.size();
    vector<vector<double> > AF(n, vector<double>(n + 1));

    for (int i = 0; i < n; ++i) {
        std::copy(A[i].begin(), A[i].end(), AF[i].begin());
        AF[i][n] = F[i];
    }

    for (int j = 0; j < n; ++j) {
        double diag = AF[j][j];
        for (int k = 0; k <= n; ++k)
            AF[j][k] /= diag;
        for (int i = 0; i < n; ++i) {
            if (i == j) continue;
            double factor = AF[i][j];
            for (int k = 0; k <= n; ++k)
                AF[i][k] -= factor * AF[j][k];
        }
    }

    vector<double> X(n);
    for (int i = 0; i < n; ++i)
        X[i] = AF[i][n];

    return X;
}

vector<vector<double> > Gauss2(vector<vector<double> > A) {
    const int n = A.size();
    vector<vector<double> > E(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i)
        E[i][i] = 1.0;

    for (int i = 0; i < n; ++i) {
        double diag = A[i][i];
        for (int j = 0; j < n; ++j) {
            A[i][j] /= diag;
            E[i][j] /= diag;
        }
        for (int j = i + 1; j < n; ++j) {
            double factor = A[j][i];
            for (int k = 0; k < n; ++k) {
                A[j][k] -= A[i][k] * factor;
                E[j][k] -= E[i][k] * factor;
            }
        }
    }

    for (int i = n - 1; i > 0; --i) {
        for (int j = i - 1; j >= 0; --j) {
            double factor = A[j][i];
            for (int k = 0; k < n; ++k) {
                A[j][k] -= A[i][k] * factor;
                E[j][k] -= E[i][k] * factor;
            }
        }
    }

    return E;
}

bool LUdecompose(const vector<vector<double> > &A, vector<vector<double> > &LU) {
    const int n = A.size();
    LU = A;

    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < i; ++k)
                sum += LU[i][k] * LU[k][j];
            LU[i][j] = A[i][j] - sum;
        }

        for (int j = i + 1; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < i; ++k)
                sum += LU[j][k] * LU[k][i];
            if (fabs(LU[i][i]) < 1e-12) return false;
            LU[j][i] = (A[j][i] - sum) / LU[i][i];
        }
    }

    return true;
}

void printLfromLU(const vector<vector<double> > &LU) {
    const int n = LU.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (j < i)
                cout << setw(15) << formatNumber(LU[i][j]);
            else if (j == i)
                cout << setw(15) << "1.00000";
            else
                cout << setw(15) << "0.00000";
        }
        cout << endl;
    }
}

void printUfromLU(const vector<vector<double> > &LU) {
    const int n = LU.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (j >= i)
                cout << setw(15) << formatNumber(LU[i][j]);
            else
                cout << setw(15) << "0.00000";
        }
        cout << endl;
    }
}

void checkLUDecomposition(const vector<vector<double> > &A, const string &label) {
    vector<vector<double> > LU;
    if (!LUdecompose(A, LU)) {
        cout << "LU-разложение " << label << " невозможно (нулевой элемент на диагонали без pivoting).\n";
        return;
    }

    cout << "LU-разложение " << label << ":\n";
    cout << "Матрица L (из " << label << "):\n";
    printLfromLU(LU);
    cout << "Матрица U (из " << label << "):\n";
    printUfromLU(LU);

    const int n = LU.size();
    vector<vector<double> > L(n, vector<double>(n, 0.0));
    vector<vector<double> > U(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            if (i > j) L[i][j] = LU[i][j];
            else if (i == j) {
                L[i][j] = 1.0;
                U[i][j] = LU[i][j];
            } else U[i][j] = LU[i][j];
        }

    auto A_reconstructed = matMult(L, U);
    cout << "Проверка: L * U (восстановленная " << label << "):\n";
    printMatrix(A_reconstructed, 1e-10);
}

vector<double> forwardSubstitution(const vector<vector<double> > &L, const vector<double> &b) {
    const int n = L.size();
    vector<double> y(n, 0.0);

    for (int i = 0; i < n; ++i) {
        y[i] = b[i];
        for (int j = 0; j < i; ++j) {
            y[i] -= L[i][j] * y[j];
        }
        y[i] /= L[i][i];
    }

    return y;
}

vector<double> backwardSubstitution(const vector<vector<double> > &U, const vector<double> &y) {
    const int n = U.size();
    vector<double> x(n, 0.0);

    for (int i = n - 1; i >= 0; --i) {
        x[i] = y[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }

    return x;
}

vector<double> LUsolve(const vector<vector<double> > &A, const vector<double> &b) {
    vector<vector<double> > LU = A;
    vector<vector<double> > L(A.size(), vector<double>(A.size(), 0.0));
    vector<vector<double> > U(A.size(), vector<double>(A.size(), 0.0));

    if (!LUdecompose(A, LU)) {
        throw std::runtime_error("LU-разложение невозможно (нулевой элемент на диагонали без pivoting).");
    }

    const int n = LU.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i > j) L[i][j] = LU[i][j];
            else if (i == j) {
                L[i][j] = 1.0;
                U[i][j] = LU[i][j];
            } else U[i][j] = LU[i][j];
        }
    }

    vector<double> y = forwardSubstitution(L, b);

    vector<double> x = backwardSubstitution(U, y);

    return x;
}


int main() {
    cout << std::scientific << std::setprecision(5);
#ifdef _WIN32
    system("chcp 65001");
#endif

    // --- Первая СЛАУ ---
    vector<vector<double> > A1 = {
        {177.20, 1.05, -8.97, 0.75},
        {4.26, 185.80, 0.13, -8.86},
        {-3.81, 5.23, -189.00, -4.88},
        {5.82, 3.87, -2.47, 81.40}
    };
    vector<double> F1 = {455.34, -924.04, -1554.46, 59.75};
    vector<double> D1 = {3, -5, 8, 1};

    // --- Вторая СЛАУ ---
    vector<vector<double> > A2 = {
        {28.859, -0.008, 2.406, 19.240},
        {14.436, -0.001, 1.203, 9.624},
        {120.204, -0.032, 10.024, 80.144},
        {-57.714, 0.016, -4.812, -38.478}
    };
    vector<double> F2 = {30.459, 18.248, 128.156, -60.908};
    vector<double> D2 = {1, 1000, -20, 3};

    // Решение
    cout << "\nСЛАУ 1 - Решение методом Гаусса:" << endl;
    vector<double> X1 = Gauss(A1, F1);
    printVector(X1);

    cout << "\nСЛАУ 2 - Решение методом Гаусса:" << endl;
    vector<double> X2 = Gauss(A2, F2);
    printVector(X2);

    // LU-разложение
    checkLUDecomposition(A1, "A1");
    checkLUDecomposition(A2, "A2");

    // LU-решения
    vector<double> x1 = LUsolve(A1, F1);
    cout << "Решение системы для A1: \n";
    printVector(x1);
    vector<double> x2 = LUsolve(A2, F2);
    cout << "Решение системы для A2: \n";
    printVector(x2);

    // Невязки
    auto R1 = vectorDiff(F1, matVecMult(A1, X1));
    auto R2 = vectorDiff(F2, matVecMult(A2, X2));

    cout << "\nНевязка для матрицы 1:" << endl;
    printVector(R1);

    cout << "\nНевязка для матрицы 2:" << endl;
    printVector(R2);

    // Нормы невязки
    cout << "\nНормы невязки СЛАУ 1:" << endl;
    cout << "  Единичная норма: " << norm1_vector(R1) << endl;
    cout << "  Бесконечная норма: " << normInf_vector(R1) << endl;

    cout << "\nНормы невязки СЛАУ 2:" << endl;
    cout << "  Единичная норма: " << norm1_vector(R2) << endl;
    cout << "  Бесконечная норма: " << normInf_vector(R2) << endl;

    // Погрешности
    auto diff1 = vectorDiff(D1, X1);
    auto diff2 = vectorDiff(D2, X2);

    cout << "\nПогрешности для СЛАУ 1:" << endl;
    double abs1 = norm1_vector(diff1), absInf1 = normInf_vector(diff1);
    cout << "  Абсолютная по 1-норме: " << abs1 << endl;
    cout << "  Абсолютная по беск. норме: " << absInf1 << endl;
    cout << "  Относительная по 1-норме: " << abs1 / norm1_vector(D1) << endl;
    cout << "  Относительная по беск. норме: " << absInf1 / normInf_vector(D1) << endl;

    cout << "\nПогрешности для СЛАУ 2:" << endl;
    double abs2 = norm1_vector(diff2), absInf2 = normInf_vector(diff2);
    cout << "  Абсолютная по 1-норме: " << abs2 << endl;
    cout << "  Абсолютная по беск. норме: " << absInf2 << endl;
    cout << "  Относительная по 1-норме: " << abs2 / norm1_vector(D2) << endl;
    cout << "  Относительная по беск. норме: " << absInf2 / normInf_vector(D2) << endl;

    // Обратные матрицы
    cout << "\nОбратная матрица A1^-1:" << endl;
    auto A1_inv = Gauss2(A1);
    printMatrix(A1_inv);
    cout << "\nПроверка A1^-1 * A1:" << endl;
    printMatrix(matMult(A1, A1_inv));

    cout << "\nОбратная матрица A2^-1:" << endl;
    auto A2_inv = Gauss2(A2);
    printMatrix(A2_inv);
    cout << "\nПроверка A2^-1 * A2:" << endl;
    printMatrix(matMult(A2, A2_inv), 1e-8);

    // Числа обусловленности
    double cond1_1 = norm1_matrix(A1) * norm1_matrix(A1_inv);
    double cond1_inf = normInf_matrix(A1) * normInf_matrix(A1_inv);

    double cond2_1 = norm1_matrix(A2) * norm1_matrix(A2_inv);
    double cond2_inf = normInf_matrix(A2) * normInf_matrix(A2_inv);

    cout << "\nЧисла обусловленности СЛАУ 1:" << endl;
    cout << "  По 1-норме: " << formatNumber(cond1_1) << endl;
    cout << "  По беск. норме: " << formatNumber(cond1_inf) << endl;

    cout << "\nЧисла обусловленности СЛАУ 2:" << endl;
    cout << "  По 1-норме: " << formatNumber(cond2_1) << endl;
    cout << "  По беск. норме: " << formatNumber(cond2_inf) << endl;

    return 0;
}
