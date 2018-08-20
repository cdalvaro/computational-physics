//
//  lu.hpp
//  NumericalPDEs
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
                
                template <class MatrixT, class MatrixC = MatrixT>
                class LU {
                public:
                    
                    typedef typename std::enable_if<std::is_floating_point<typename MatrixT::ValueType>::value >::type *IsFloatingPoint;
                    typedef typename MatrixT::ValueType ValueType;
                    
                    LU(const MatrixC &matrix) :
                    is_factorized(false), original(matrix), rows(matrix.Rows()) {
                        if (!matrix.IsSquared()) {
                            throw std::logic_error("LU matrix cannot be computed for a non-squared matrix.");
                        }
                    }
                    
                    virtual ~LU() = default;
                    
                    const MatrixT &L() {
                        FactorizeLU();
                        return l;
                    }
                    
                    const MatrixT &U() {
                        FactorizeLU();
                        return u;
                    }
                    
                    template <class VectorC>
                    VectorC SolveLinearSystem(const VectorC &b_terms) {
                        if (rows != b_terms.Size()) {
                            throw std::logic_error("The number of rows of the LU matrix does not match the number of elements in the b terms vector.");
                        }
                        
                        FactorizeLU();
                        
                        VectorC tmp(rows, 0), x(rows, 0);
                        
                        ValueType sum;
                        for (size_t row = 0; row < rows; ++row) {
                            sum = 0;
                            for (size_t column = 0; column < row; ++column) {
                                sum += original[row][column] * tmp[column];
                            }
                            tmp[row] = b_terms[row] - sum;
                        }
                        
                        // Since matrix is square -> #rows == #columns
                        const size_t last_row = rows - 1;
                        for (ssize_t row = last_row; row >= 0; --row) {
                            sum = 0;
                            for (ssize_t column = last_row; column >= 0; --column) {
                                sum += original[row][column] * x[column];
                            }
                            x[row] = (tmp[row] - sum) / original[row][row];
                        }
                        
                        return x;
                    }
                    
                    MatrixC InverseMatrix(IsFloatingPoint = nullptr) {
                        
                        MatrixC inverse(rows, rows);
                        const auto I = MatrixC::Identity(rows);
                        
                        for (size_t k = 0; k < rows; ++k) {
                            inverse.SetColumn(k, SolveLinearSystem(I.GetRowAsVector(k)));
                        }
                        
                        return inverse;
                    }
                    
                    static MatrixC InverseMatrix(const MatrixC &matrix, IsFloatingPoint = nullptr) {
                        return LU<MatrixT, MatrixC>(matrix).InverseMatrix();
                    }
                    
                    ValueType Determinant(IsFloatingPoint = nullptr) {
                        FactorizeLU();
                        
                        ValueType determinant = u[0][0];
                        for (size_t row = 1; row < rows; ++row) {
                            determinant *= u[row][row];
                        }
                        
                        return determinant;
                    }
                    
                    static ValueType Determinant(const MatrixC &matrix, IsFloatingPoint = nullptr) {
                         return LU<MatrixT, MatrixC>(matrix).Determinant();
                    }
                    
                private:
                    
                    const MatrixC original;
                    const size_t rows;
                    
                    MatrixT l, u;
                    
                    bool is_factorized;
                    
                    void FactorizeLU() {
                        
                        if (is_factorized) {
                            return;
                        }
                        
                        MatrixT lu(original);
                        
                        auto first_element = lu[0][0];
                        if (first_element == 0) {
                            throw std::logic_error("Unable to factorize LU matrix because first element is 0");
                        }
                        
                        for (size_t row = 1; row < rows; ++row) {
                            lu[row][0] /= first_element;
                        }

                        ValueType sum = 0;
                        const size_t last_row = rows - 1;
                        for (size_t row = 1; row < last_row; ++row) {

                            // U Matrix
                            for (size_t column = row; column < rows; ++column) {
                                sum = 0;
                                for (size_t k = 0; k < row; ++k) {
                                    sum += lu[k][column] * lu[row][k];
                                }
                                lu[row][column] -= sum;
                            }

                            // L Matrix
                            for (size_t k = row + 1; k < rows; ++k) {
                                sum = 0;
                                for (size_t column = 0; column < row; ++column) {
                                    sum += lu[column][row] * lu[k][column];
                                }
                                lu[k][row] = (lu[k][row] - sum) / lu[row][row];
                            }
                        }

                        sum = 0;
                        for (size_t k = 0; k < last_row; ++k) {
                            sum += lu[k][last_row] * lu[last_row][k];
                        }
                        lu[last_row][last_row] -= sum;

                        l.Resize(rows, rows, 0);
                        u.Resize(rows, rows, 0);
                        
                        for (size_t row = 0; row < rows; ++row) {
                            for (size_t column = 0; column < rows; ++column) {
                                if (column > row) {
                                    u[row][column] = lu[row][column];
                                    l[row][column] = 0;
                                } else if (column == row) {
                                    u[row][column] = lu[row][column];
                                    l[row][column] = 1;
                                } else {
                                    u[row][column] = 0;
                                    l[row][column] = lu[row][column];
                                }
                            }
                        }
                        
                        is_factorized = true;
                    }
                    
                };
                
            } /* namespace factorization */
        } /* namespace algorithms */
    } /* namespace math */
} /* namespace cda */
