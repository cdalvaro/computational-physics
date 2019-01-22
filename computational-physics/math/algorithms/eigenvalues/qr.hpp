//
//  qr.hpp
//  Computational Physics
//
//  Created by Carlos David on 17/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#pragma once

#include <map>

#include "../factorization/lu.hpp"
#include "../../containers/vector.hpp"
#include "../../math.hpp"


#define CDA_QR_DEFAULT_ACCURACY 1E-06
#define CDA_QR_DEFAULT_MAX_ITERATIONS 1E+05


namespace cda {
    namespace math {
        namespace algorithms {
            namespace eigenvalues {
                
                template <template<typename T> class Matrix, typename ValueType = double,
                          class = typename std::enable_if<std::is_floating_point<ValueType>::value>::type>
                class QR {
                public:
                    
                    QR(const Matrix<ValueType> &matrix,
                       const double &accuracy = CDA_QR_DEFAULT_ACCURACY,
                       const size_t &max_iterations = CDA_QR_DEFAULT_MAX_ITERATIONS) :
                    original(matrix), rows(matrix.Rows()),
                    max_iterations(max_iterations), accuracy(accuracy) {
                        if (!matrix.IsSquare()) {
                            throw std::logic_error("Matrix must be square to compute its eigenvalues.");
                        }
                    }
                    
                    virtual ~QR() = default;
                    
                    const Matrix<ValueType> &Q() {
                        if (q.IsNull()) {
                            ComputeQR(original);
                        }
                        return q;
                    }
                    
                    const Matrix<ValueType> &R() {
                        if (r.IsNull()) {
                            ComputeQR(original);
                        }
                        return r;
                    }
                    
                    const size_t &MaxIterations() const {
                        return max_iterations;
                    }
                    
                    void MaxIterations(const size_t &max_iterations) {
                        this->max_iterations = max_iterations;
                    }
                    
                    const double &Accuracy() const {
                        return accuracy;
                    }
                    
                    void Accuracy(const double &accuracy) {
                        this->accuracy = accuracy;
                    }
                    
                    const containers::Vector<ValueType> &EigenValues() {
                        if (eigen_values.IsEmpty()) {
                            auto matrix(original);
                            const size_t last_row = rows - 1;
                            ValueType square_sum, element;
                            
                            for (size_t iteration = 0; iteration < max_iterations; ++iteration) {
                                ComputeQR(matrix);
                                matrix = r * q;
                                
                                // Convergence test
                                square_sum = 0;
                                for (size_t column = 0; column < last_row; ++column) {
                                    for (size_t row = column + 1; row < rows; ++row) {
                                        element = matrix[row][column];
                                        square_sum += element * element;
                                    }
                                }
                                
                                if (std::sqrt(square_sum) < accuracy) {
                                    break;
                                }
                            }
                            
                            eigen_values = matrix.GetDiagonal();
                        }
                        
                        return eigen_values;
                    }
                    
                    const containers::Vector<ValueType> &EigenVector(const ValueType &eigen_value) {
                        
                        auto it_eigen_vector = eigen_vectors.find(eigen_value);
                        if (it_eigen_vector != eigen_vectors.end() && !it_eigen_vector->second.IsEmpty()) {
                            return it_eigen_vector->second;
                        }
                        
                        Matrix<ValueType> inverse_matrix(rows, rows, 0);
                        inverse_matrix.SetDiagonal(eigen_value * (accuracy + 1.0));
                        
                        inverse_matrix = (original - inverse_matrix).Pow(-1);
                        
                        ValueType normalization_factor = 0.0;
                        ValueType old_normalization_factor, distance;
                        Matrix<ValueType> eigenVector(rows, 1, 1);
                        
                        for (size_t iteration = 0; iteration < max_iterations; ++iteration) {
                            old_normalization_factor = normalization_factor;
                            
                            eigenVector = inverse_matrix * eigenVector;
                            normalization_factor = eigenVector.abs_max_element_with_sign();
                            eigenVector /= normalization_factor;
                            
                            // Convergence test
                            distance = normalization_factor - old_normalization_factor;
                            if (std::sqrt(distance * distance) < accuracy) {
                                break;
                            }
                        }
                        
                        eigen_vectors.emplace(eigen_value, eigenVector.GetColumnAsVector(0) / eigenVector[0][rows - 1]);
                        
                        return eigen_vectors[eigen_value];
                    }
                    
                    const std::map<ValueType, containers::Vector<ValueType>> &EigenVectors() {
                        auto values = EigenValues();
                        for (auto it_value = values.begin(); it_value != values.end(); ++it_value) {
                            EigenVector(*it_value);
                        }
                        return eigen_vectors;
                    }
                    
                private:
                    
                    const Matrix<ValueType> original;
                    const size_t rows;
                    
                    size_t max_iterations;
                    double accuracy;
                    
                    Matrix<ValueType> q, r;
                    containers::Vector<ValueType> eigen_values;
                    std::map<ValueType, containers::Vector<ValueType>> eigen_vectors;
                    
                    void ComputeQR(const Matrix<ValueType> &matrix) {
                        
                        const auto I = Matrix<ValueType>::Identity(rows);
                        
                        //  First column
                        auto c = matrix.GetColumnAsVector(0);
                        auto vt = I.GetRowAsVector(0);
                        vt *= signum(c[0]) * c.Norm();
                        vt += c;
                        
                        q = containers::Transpose(vt) * vt;
                        q *= -2.0 / vt.SquareNorm();
                        q += I;
                        
                        r = q * matrix;
                        for (auto it_r = r.begin() + rows; it_r != r.end(); it_r += rows) {
                            *it_r = 0;
                        }
                        
                        //  Other columns
                        Matrix<ValueType> h, h__;
                        const auto last_row = rows - 1;
                        ValueType *it_r_column, *it_end_r_column;
                        
                        for (size_t row = 1; row < last_row; ++row) {
                            c = r.GetColumnAsVector(row, row);
                            
                            if (c.IsNull()) {
                                h = I;
                            } else {
                                vt = I.GetRowAsVector(row, row);
                                vt *= signum(c[0]) * c.Norm();
                                vt += c;
                                
                                h = containers::Transpose(vt) * vt;
                                h *= -2.0 / vt.SquareNorm();
                                h += Matrix<ValueType>::Identity(rows - row);
                                
                                h__ = I;
                                h__.SetMatrix(row, row, h);
                                h = std::move(h__);
                            }
                            
                            q = h * q;
                            r = h * r;
                            
                            it_r_column = r[row];
                            it_end_r_column = it_r_column + row;
                            for ( ; it_r_column != it_end_r_column; ++it_r_column) {
                                *it_r_column = 0;
                            }
                        }
                        
                        q = std::move(q.Transpose());
                        
                        for (size_t row = 1; row < rows; ++row) {
                            it_r_column = r[row];
                            it_end_r_column = it_r_column + row;
                            for ( ; it_r_column < it_end_r_column; ++it_r_column) {
                                *it_r_column = 0;
                            }
                        }
                        
                    }
                    
                };
                
            } /* namespace eigenvalues */
        } /* namespace algorithms */
    } /* namespace math */
} /* namespace cda */
