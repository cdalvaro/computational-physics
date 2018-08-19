//
//  lu.hpp
//  NumericalPDEs
//
//  Created by Carlos David on 16/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#pragma once

#include "../../containers.hpp"


namespace cda {
    namespace math {
        namespace algorithms {
            namespace factorization {
                
                template <typename T>
                class LU {
                public:
                    
                    LU(const containers::Matrix<T> &matrix) :
                    is_factorized(false), original(matrix), rows(matrix.Rows()) {
                        if (!matrix.IsSquared()) {
                            std::logic_error("LU matrix cannot be computed for a non-squared matrix.");
                        }
                    }
                    
                    virtual ~LU() = default;
                    
                    const containers::Matrix<T> &L() {
                        FactorizeLU();
                        return l;
                    }
                    
                    const containers::Matrix<T> &U() {
                        FactorizeLU();
                        return u;
                    }
                    
                    containers::Vector<T> SolveLinearSystem(const containers::Vector<T> &b_terms) {
                        if (rows != b_terms.Size()) {
                            throw std::logic_error("The number of rows of the LU matrix does not match the number of elements in the b terms vector.");
                        }
                        
                        FactorizeLU();
                        
                        containers::Vector<T> tmp(rows, 0), x(rows, 0);
                        
                        T sum;
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
                    
                    containers::Matrix<T> InverseMatrix() {
                        
                        containers::Matrix<T> inverse(rows, rows);
                        containers::Vector<T> I(rows, 0);
                        
                        for (size_t k = 0; k < rows; ++k) {
                            I[k] = 1;
                            inverse.SetColumn(k, SolveLinearSystem(I));
                            I[k] = 0;
                        }
                        
                        return inverse;
                    }
                    
                    static containers::Matrix<T> InverseMatrix(const containers::Matrix<T> &matrix) {
                        return LU<T>(matrix).InverseMatrix();
                    }
                    
                    T Determinant() const {
                        T determinant = 1.0;
                        
                        FactorizeLU();
                        for (size_t row = 0; row < rows; ++row) {
                            determinant *= u[row][row];
                        }
                        
                        return determinant;
                    }
                    
                    static T Determinant(const containers::Matrix<T> &matrix) {
                         return LU<T>(matrix).Determinant();
                    }
                    
                private:
                    
                    const containers::Matrix<T> original;
                    const size_t rows;
                    
                    containers::Matrix<T> l;
                    containers::Matrix<T> u;
                    
                    bool is_factorized;
                    
                    void FactorizeLU() {
                        
                        if (is_factorized) {
                            return;
                        }
                        
                        auto matrix(original);
                        
                        l.Resize(rows, rows, 0);
                        u.Resize(rows, rows, 0);
                        
                        auto first_element = matrix[0][0];
                        for (size_t row = 1; row < rows; ++row) {
                            matrix[row][0] /= first_element;
                        }
                        
                        T sum;
                        const size_t last_row = rows - 1;
                        for (size_t row = 1; row < last_row; ++row) {
                            
                            // U Matrix
                            for (size_t column = row; column < rows; ++column) {
                                sum = 0;
                                for (size_t k = 0; k < row; ++k) {
                                    sum += matrix[k][column] * matrix[row][k];
                                }
                                matrix[row][column] -= sum;
                            }
                            
                            // L Matrix
                            for (size_t k = row + 1; k < rows; ++k) {
                                sum = 0;
                                for (size_t column = 0; column < row; ++column) {
                                    sum += matrix[column][row] * matrix[k][column];
                                }
                                matrix[k][row] = (matrix[k][row] - sum) / matrix[row][row];
                            }
                        }
                        
                        sum = 0;
                        for (size_t k = 0; k < last_row; ++k) {
                            sum += matrix[k][last_row] * matrix[last_row][k];
                        }
                        matrix[last_row][last_row] -= sum;
                        
                        for (size_t row = 0; row < rows; ++row) {
                            for (size_t column = 0; column < rows; ++column) {
                                if (column > row) {
                                    u[row][column] = matrix[row][column];
                                    l[row][column] = 0;
                                } else if (column == row) {
                                    u[row][column] = matrix[row][column];
                                    l[row][column] = 1;
                                } else {
                                    u[row][column] = 0;
                                    l[row][column] = matrix[row][column];
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
