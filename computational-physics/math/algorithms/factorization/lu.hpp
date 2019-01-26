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
                    lu(matrix), rows(matrix.rows()), is_factorized(false), _is_degenerate(false) {
                        if (!matrix.is_square()) {
                            throw std::logic_error("LU matrix cannot be computed for a non-square matrix.");
                        }
                    }
                    
                    virtual ~LU() = default;
                    
                    const Matrix<ValueType> &l() {
                        FactorizeLU();
                        return _l;
                    }
                    
                    const Matrix<ValueType> &u() {
                        FactorizeLU();
                        return _u;
                    }
                    
                    const bool &is_degenerate() const {
                        return _is_degenerate;
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
                        if (_is_degenerate) {
                            throw std::logic_error("Matrix is degenerate, so does not have inverse");
                        }
                        
                        Matrix<ValueType> inverse(rows, rows);
                        const auto I = Matrix<ValueType>::identity(rows);
                        
                        for (size_t k = 0; k < rows; ++k) {
                            inverse.set_column(k, SolveLinearSystem(I.get_row_as_vector(k)));
                        }
                        
                        return inverse;
                    }
                    
                    template<typename OtherType>
                    static Matrix<ValueType> InverseMatrix(const Matrix<OtherType> &matrix) {
                        return LU<Matrix, ValueType>(matrix).InverseMatrix();
                    }
                    
                    ValueType determinant() {
                        
                        FactorizeLU();
                        if (_is_degenerate) {
                            return 0;
                        }
                        
                        ValueType determinant = _u[0][0];
                        for (size_t row = 1; row < rows; ++row) {
                            determinant *= _u[row][row];
                        }
                        
                        return determinant;
                    }
                    
                    template<typename OtherType>
                    static ValueType determinant(const Matrix<OtherType> &matrix) {
                         return LU<Matrix, ValueType>(matrix).determinant();
                    }
                    
                private:
                    
                    Matrix<ValueType> _l, _u, lu;
                    const size_t rows;
                    
                    bool is_factorized;
                    bool _is_degenerate;
                    
                    void FactorizeLU() {
                        if (is_factorized) {
                            return;
                        }
                        
                        RemoveSingularities();
                        
                        _l.resize(rows, rows, 0);
                        _u.resize(rows, rows, 0);
                        
                        for (size_t row = 0; row < rows; ++row) {
                            _l[row][row] = 1;
                        }
                        
                        for (size_t row = 0; row < rows; ++row) {
                            for (size_t column = 0; column < rows; ++column) {
                                double sum;
                                if (column <= row) {
                                    sum = 0;
                                    for (size_t k = 0; k < column; ++k) {
                                        sum += _l[column][k] * _u[k][row];
                                    }
                                    _u[column][row] = lu[column][row] - sum;
                                }
                                
                                if (column >= row) {
                                    sum = 0;
                                    for (size_t k = 0; k < row; ++k) {
                                        sum += _l[column][k] * _u[k][row];
                                    }
                                    _l[column][row] = (lu[column][row] - sum) / _u[row][row];
                                }
                            }
                        }
                        
                        for (size_t row = 0; row < rows; ++row) {
                            for (size_t column = 0; column < rows; ++column) {
                                if (column > row) {
                                    lu[row][column] = _u[row][column];
                                    _l[row][column]  = 0;
                                } else if (column == row) {
                                    lu[row][column] = _u[row][column];
                                    _l[row][column]  = 1;
                                    
                                    if (!_is_degenerate && _u[row][column] == 0) {
                                        _is_degenerate = true;
                                    }
                                } else {
                                    _u[row][column]  = 0;
                                    lu[row][column] = _l[row][column];
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
