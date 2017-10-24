//
// Created by gritukan on 10/24/17.
//

#ifndef GAUSS_MATRIX_H
#define GAUSS_MATRIX_H

#include <vector>
#include <functional>

template<typename T>
class Matrix {
private:
    std::vector<std::vector<T> > matrix;

    static T Tabs(T x) {
        return (x < T::ZERO) ? negate(x) : x;
    }

public:
    Matrix() = default;
    ~Matrix() = default;

    Matrix(size_t n, size_t m) {
        matrix.resize(n, std::vector<T>(m));
    }

    Matrix(const Matrix &other) = default;
    Matrix& operator=(const Matrix &other) = default;
    Matrix(Matrix &&other) noexcept = default;

    size_t getN() const {
        return matrix.size();
    }

    size_t getM() const {
        return matrix[0].size();
    }

    decltype(matrix.at(0)) operator[](size_t n) {
        return matrix.at(n);
    }

    const decltype(matrix.at(0)) operator[](size_t n) const {
        return matrix.at(n);
    }

    template<typename U>
    friend void makeGaussianElimination(Matrix<U> &matrix, bool selectSmallest,
                                        std::function<void(size_t, size_t)> onSwap,
                                        std::function<void(size_t, U)> onDivide,
                                        std::function<void(size_t, size_t, U)> onSubtract);

    T determant(bool selectSmallest) const {
        if (getN() != getM()) {
            throw std::logic_error("Matrix should have square form");
        }
        Matrix toSolve(*this);
        bool swapsParity = false;

        std::function<void(size_t, size_t)> onSwap = [&] (size_t x, size_t y) -> void {
            if (x != y) swapsParity ^= true;
        };
        T result = T::ONE;
        std::function<void(size_t, T)> onDivide = [&] (size_t x, T value) -> void {
            result = result * value;
        };
        std::function<void(size_t, size_t, T)> onSubstract = [&] (size_t x, size_t y, T value) -> void {};

        makeGaussianElimination(toSolve, selectSmallest, onSwap, onDivide, onSubstract);
        for (size_t i = 0; i < getN(); i++) {
            result = result * toSolve.matrix[i][i];
        }
        if (swapsParity) {
            result = negate(result);
        }
        return result;
    };


};

template<typename T>
void makeGaussianElimination(Matrix<T> &matrix, bool selectSmallest,
                              std::function<void(size_t, size_t)> onSwap,
                              std::function<void(size_t, T)> onDivide,
                              std::function<void(size_t, size_t, T)> onSubtract) {
    size_t n = matrix.getN(), m = matrix.getM();
    auto &a = matrix.matrix;
    size_t nextColumn = 0;

    for (size_t i = 0; i < n; i++) {
        size_t position = n;
        while (nextColumn < m && position == n) {
            for (size_t j = i; j < n; j++) {
                if (!(a[j][nextColumn] == T::ZERO) &&
                        (!selectSmallest || position == n ||
                         Matrix<T>::Tabs(a[j][nextColumn]) > Matrix<T>::Tabs(a[position][nextColumn]))) {
                    position = j;
                }
            }
            nextColumn++;
        }
        if (position != n) {
            size_t currentColumn = nextColumn - 1;
            onSwap(position, i);
            if (position != i) {
                std::swap(a[i], a[position]);
            }
            onDivide(i, a[i][currentColumn]);
            T divisor = a[i][currentColumn];
            for (size_t j = currentColumn; j < m; j++) {
                a[i][j] = a[i][j] * inverse(divisor);
            }
            for (size_t j = i + 1; j < n; j++) {
                T multiplier = a[j][currentColumn];
                onSubtract(j, i, multiplier);
                for (size_t k = currentColumn; k < m; k++) {
                    a[j][k] = a[j][k] + negate(a[i][k] * multiplier);
                }
            }
        }
    }
}


#endif //GAUSS_MATRIX_H
