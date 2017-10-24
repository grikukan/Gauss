#include <iostream>
#include "FieldNumber.h"
#include "RealNumber.h"
#include "Matrix.h"
int main() {
    int n;
    std::cin >> n;
    Matrix<RealNumber<long double> > m(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double x;
            std::cin >> x;
            m[i][j] = x;
        }
    }
    std::cout.precision(20);

    long double res = m.determant(false);
    std::cout << std::fixed << res << std::endl;
}