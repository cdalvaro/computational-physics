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
                    is_factorized(false) {
                        
                        if (!matrix.IsSquared()) {
                            std::logic_error("LU matrix cannot be computed for a non-squared matrix.");
                        }
                        
                        lu = matrix;
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
                        
                        const auto rows = lu.Rows();
                        const auto columns = lu.Columns();
                        
                        if (rows != b_terms.Size()) {
                            throw std::logic_error("The number of rows of the LU matrix does not match the number of elements in the b terms vector.");
                        }
                        
                        FactorizeLU();
                        
                        containers::Vector<T> tmp(rows, 0), x(rows, 0);
                        
                        T sum;
                        for (size_t row = 0; row < rows; ++row) {
                            sum = 0;
                            for (size_t column = 0; column < row; ++column) {
                                sum += lu[row][column] * tmp[column];
                            }
                            tmp[row] = b_terms[row] - sum;
                        }
                        
                        for (ssize_t row = rows - 1; row >= 0; --row) {
                            sum = 0;
                            for (ssize_t column = columns - 1; column >= 0; --column) {
                                sum += lu[row][column] * x[column];
                            }
                            x[row] = (tmp[row] - sum) / lu[row][row];
                        }
                        
                        return x;
                    }
                    
                    containers::Matrix<T> InverseMatrix() {
                        
                        const auto rows = lu.Rows();
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
                        LU<T> lu_matrix(matrix);
                        return lu_matrix.InverseMatrix();
                    }
                    
                private:
                    
                    containers::Matrix<T> lu;
                    containers::Matrix<T> l;
                    containers::Matrix<T> u;
                    
                    bool is_factorized;
                    
                    void FactorizeLU() {
                        
                        if (is_factorized) {
                            return;
                        }
                        
                        auto rows = lu.Rows();
                        auto columns = lu.Columns();
                        
                        l.Resize(rows, columns, 0);
                        u.Resize(rows, columns, 0);
                        
                        auto first_element = lu[0][0];
                        for (size_t row = 1; row < rows; ++row) {
                            lu[row][0] /= first_element;
                        }
                        
                        T sum;
                        for (size_t row = 1; row < rows - 1; ++row) {
                            
                            // U Matrix
                            for (size_t column = row; column < columns; ++column) {
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
                        for (size_t k = 0; k < rows - 1; ++k) {
                            sum += lu[k][columns - 1] * lu[rows - 1][k];
                        }
                        lu[rows - 1][columns - 1] -= sum;
                        
                        for (size_t row = 0; row < rows; ++row) {
                            for (size_t column = 0; column < columns; ++column) {
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
