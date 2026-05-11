#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>

using namespace std;

#define MAX_SIZE 500

// Tính autocorrelation của signal
void computeAutocorrelation(double signal[MAX_SIZE], double autocorr[MAX_SIZE], int N) {
    // R_xx(k) = (1/N) * sum_{n=k}^{N-1} x[n] * x[n-k]
    for (int k = 0; k < N; k++) {
        double sum = 0.0;
        for (int n = k; n < N; n++)
            sum += signal[n] * signal[n - k];
        autocorr[k] = sum / N;
    }
}

// Hàm tính crosscorrelation giữa desired signal và input signal
void computeCrosscorrelation(double desired[MAX_SIZE], double input[MAX_SIZE], double crosscorr[MAX_SIZE], int N) {
    // p_dx(l) = (1/N) * sum_{n=l}^{N-1} d[n] * x[n-l]
    for (int l = 0; l < N; l++) {
        double sum = 0.0;
        for (int n = l; n < N; n++)
            sum += desired[n] * input[n - l];
        crosscorr[l] = sum / N;
    }
}

// Hàm tạo ma trận Toeplitz từ autocorrelation
void createToeplitzMatrix(double autocorr[MAX_SIZE], double R[MAX_SIZE][MAX_SIZE], int N) {
    // Symmetric Toeplitz: R[i][j] = R_xx(|i - j|)
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            R[i][j] = autocorr[abs(i - j)];
}

// Giải hệ phương trình tuyến tính bằng Gauss elimination
void solveLinearSystem(double A[MAX_SIZE][MAX_SIZE], double b[MAX_SIZE], double x[MAX_SIZE], int N) {
    // Gaussian elimination with partial pivoting
    // Use static to avoid stack overflow (500x500 doubles = 2MB)
    static double a[MAX_SIZE][MAX_SIZE];
    double rhs[MAX_SIZE];
    for (int i = 0; i < N; i++) {
        rhs[i] = b[i];
        for (int j = 0; j < N; j++)
            a[i][j] = A[i][j];
    }

    for (int col = 0; col < N; col++) {
        // Partial pivot: find row with largest |a[row][col]|
        int pivot = col;
        for (int row = col + 1; row < N; row++)
            if (fabs(a[row][col]) > fabs(a[pivot][col]))
                pivot = row;

        // Swap rows
        for (int j = 0; j < N; j++) swap(a[col][j], a[pivot][j]);
        swap(rhs[col], rhs[pivot]);

        // Eliminate below
        for (int row = col + 1; row < N; row++) {
            if (fabs(a[col][col]) < 1e-15) continue;
            double factor = a[row][col] / a[col][col];
            for (int k = col; k < N; k++)
                a[row][k] -= factor * a[col][k];
            rhs[row] -= factor * rhs[col];
        }
    }

    // Back substitution
    for (int i = N - 1; i >= 0; i--) {
        x[i] = rhs[i];
        for (int j = i + 1; j < N; j++)
            x[i] -= a[i][j] * x[j];
        if (fabs(a[i][i]) > 1e-15)
            x[i] /= a[i][i];
    }
}

// Tính hệ số Wiener
void computeWienerCoefficients(double desired[MAX_SIZE], double input[MAX_SIZE], int N, double coefficients[MAX_SIZE]) {
    double autocorr[MAX_SIZE];
    double crosscorr[MAX_SIZE];
    double R[MAX_SIZE][MAX_SIZE];

    computeAutocorrelation(input, autocorr, N);
    computeCrosscorrelation(desired, input, crosscorr, N);
    createToeplitzMatrix(autocorr, R, N);
    solveLinearSystem(R, crosscorr, coefficients, N);
}

// Áp dụng Wiener filter
void applyWienerFilter(double input[MAX_SIZE], double coefficients[MAX_SIZE], double output[MAX_SIZE], int N) {
    // y[n] = sum_{k=0}^{N-1} h[k] * x[n-k]  (causal: only use x[n-k] when n-k >= 0)
    for (int n = 0; n < N; n++) {
        double sum = 0.0;
        for (int k = 0; k < N && (n - k) >= 0; k++)
            sum += coefficients[k] * input[n - k];
        output[n] = sum;
    }
}

// Tính MMSE
double computeMMSE(double desired[MAX_SIZE], double output[MAX_SIZE], int N) {
    // MMSE = mean of squared differences between desired and filtered output
    double sum = 0.0;
    for (int n = 0; n < N; n++) {
        double err = desired[n] - output[n];
        sum += err * err;
    }
    return sum / N;
}

// Đọc file
int readSignalFromFile(const string &filename, double signal[MAX_SIZE]) {
    ifstream file(filename);
    if (!file.is_open()) throw runtime_error("Cannot open file: " + filename);

    int count = 0;
    double val;
    while (file >> val)
        signal[count++] = val;
    return count;
}

// Ghi file
void writeOutputToFile(const string &filename, double output[MAX_SIZE], int N, double mmse)
{
    ofstream file(filename);
    if (!file.is_open())
        throw runtime_error("Cannot create output file");

    file << "Filtered output: ";

    file << fixed << setprecision(1);
    for (int i = 0; i < N; i++)
    {
        double val = output[i];
        val = round(val * 10.0) / 10.0;
        if (fabs(val) < 1e-9)
            val = 0.0;
        file << val;
        if (i != N - 1)
            file << " ";
    }

    file
        << "\nMMSE: " << setprecision(1) << mmse << endl;

    file.close();
}

int main()
{
    try
    {
        double desired[MAX_SIZE] = {0}, input[MAX_SIZE] = {0}, output[MAX_SIZE] = {0}, coefficients[MAX_SIZE] = {0};

        int SIZE = readSignalFromFile("desired.txt", desired);
        int N2 = readSignalFromFile("input.txt", input);

        if (SIZE != N2)
        {
            ofstream errorFile("output.txt");
            errorFile << "Error: size not match" << endl;
            errorFile.close();
            cerr << "Error: size not match" << endl;
            return 0;
        }

        computeWienerCoefficients(desired, input, SIZE, coefficients);
        applyWienerFilter(input, coefficients, output, SIZE);
        double mmse = computeMMSE(desired, output, SIZE);
        writeOutputToFile("output.txt", output, SIZE, mmse);
        cout << "Done VO TIEN ! Check output.txt for results." << endl;
    }
    catch (const exception &e)
    {
        cerr << "Error: " << e.what() << endl;
        ofstream errorFile("output.txt");
        errorFile << e.what() << endl;
        errorFile.close();
        return 0;
    }

    return 0;
}