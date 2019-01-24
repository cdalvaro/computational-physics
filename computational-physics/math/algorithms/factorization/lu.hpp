//
//  lu.hpp
//  Computational Physics
//
//  Created by Carlos David on 16/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#pragma once

#include <stdexcept>


namespace cda {
    namespace math {
        namespace algorithms {
            namespace factorization {
                
                template <template<typename T> class Matrix, typename ValueType = double,
                          class = typename std::enable_if<std::is_floating_point<ValueType>::value>::type>
                class LU {
                public:
                    
                    LU(const Matrix<ValueType> &matrix) :
                    lu(matrix), rows(matrix.Rows()), is_factorized(false), is_degenerate(false) {
                        if (!matrix.IsSquare()) {
                            throw std::logic_error("LU matrix cannot be computed for a non-square matrix.");
                        }
                    }
                    
                    virtual ~LU() = default;
                    
                    const Matrix<ValueType> &L() {
                        FactorizeLU();
                        return l;
                    }
                    
                    const Matrix<ValueType> &U() {
                        FactorizeLU();
                        return u;
                    }
                    
                    const bool &IsDegenerate() const {
                        return is_degenerate;
                    }
                    
                    template <template<typename> class Vector>
                    Vector<ValueType> SolveLinearSystem(const Vector<ValueType> &b_terms) {
                        
                        if (rows != b_terms.size()) {
                            throw std::logic_error("The number of rows of the LU matrix does not match the number of elements in the b terms vector.");
                        }
                        
                        FactorizeLU();
                        
                        Vector<ValueType> tmp(rows, 0), x(rows, 0);
                        
                        ValueType sum;
                        for (size_t row = 0; row < rows; ++row) {
                            sum = 0;
                            for (size_t column = 0; column < row; ++column) {
                                sum += lu[row][column] * tmp[column];
                            }
                            tmp[row] = b_terms[row] - sum;
                        }
                        
                        // Since matrix is square -> #rows == #columns
                        const size_t last_row = rows - 1;
                        for (ssize_t row = last_row; row >= 0; --row) {
                            sum = 0;
                            for (ssize_t column = last_row; column >= 0; --column) {
                                sum += lu[row][column] * x[column];
                            }
                            x[row] = (tmp[row] - sum) / lu[row][row];
                        }
                        
                        return x;
                    }
                    
                    Matrix<ValueType> InverseMatrix() {
                        
                        FactorizeLU();
                        if (is_degenerate) {
                            throw std::logic_error("Matrix is degenerate, so does not have inverse");
                        }
                        
                        Matrix<ValueType> inverse(rows, rows);
                        const auto I = Matrix<ValueType>::Identity(rows);
                        
                        for (size_t k = 0; k < rows; ++k) {
                            inverse.SetColumn(k, SolveLinearSystem(I.GetRowAsVector(k)));
                        }
                        
                        return inverse;
                    }
                    
                    template<typename OtherType>
                    static Matrix<ValueType> InverseMatrix(const Matrix<OtherType> &matrix) {
                        return LU<Matrix, ValueType>(matrix).InverseMatrix();
                    }
                    
                    ValueType Determinant() {
                        
                        FactorizeLU();
                        if (is_degenerate) {
                            return 0;
                        }
                        
                        ValueType determinant = u[0][0];
                        for (size_t row = 1; row < rows; ++row) {
                            determinant *= u[row][row];
                        }
                        
                        return determinant;
                    }
                    
                    template<typename OtherType>
                    static ValueType Determinant(const Matrix<OtherType> &matrix) {
                         return LU<Matrix, ValueType>(matrix).Determinant();
                    }
                    
                private:
                    
                    Matrix<ValueType> l, u, lu;
                    const size_t rows;
                    
                    bool is_factorized;
                    bool is_degenerate;
                    
                    void FactorizeLU() {
                        if (is_factorized) {
                            return;
                        }
                        
                        RemoveSingularities();
                        
                        l.resize(rows, rows, 0);
                        u.resize(rows, rows, 0);
                        
                        for (size_t row = 0; row < rows; ++row) {
                            l[row][row] = 1;
                        }
                        
                        for (size_t row = 0; row < rows; ++row) {
                            for (size_t column = 0; column < rows; ++column) {
                                double sum;
                                if (column <= row) {
                                    sum = 0;
                                    for (size_t k = 0; k < column; ++k) {
                                        sum += l[column][k] * u[k][row];
                                    }
                                    u[column][row] = lu[column][row] - sum;
                                }
                                
                                if (column >= row) {
                                    sum = 0;
                                    for (size_t k = 0; k < row; ++k) {
                                        sum += l[column][k] * u[k][row];
                                    }
                                    l[column][row] = (lu[column][row] - sum) / u[row][row];
                                }
                            }
                        }
                        
                        for (size_t row = 0; row < rows; ++row) {
                            for (size_t column = 0; column < rows; ++column) {
                                if (column > row) {
                                    lu[row][column] = u[row][column];
                                    l[row][column]  = 0;
                                } else if (column == row) {
                                    lu[row][column] = u[row][column];
                                    l[row][column]  = 1;
                                    
                                    if (!is_degenerate && u[row][column] == 0) {
                                        is_degenerate = true;
                                    }
                                } else {
                                    u[row][column]  = 0;
                                    lu[row][column] = l[row][column];
                                }
                            }
                        }
                        
                        is_factorized = true;
                    }
                    
                    void RemoveSingularities() {
                        auto diagonal = lu.get_diagonal();
                        if (diagonal.find(0) == diagonal.end()) {
                            return;
                        }
                        
                        throw std::logic_error("Matrix is nonsingular");
                    }
                    
                };
                
            } /* namespace factorization */
        } /* namespace algorithms */
    } /* namespace math */
} /* namespace cda */
