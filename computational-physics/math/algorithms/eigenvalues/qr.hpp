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
                    original(matrix), rows(matrix.rows()),
                    max_iterations(max_iterations), accuracy(accuracy) {
                        if (!matrix.is_square()) {
                            throw std::logic_error("Matrix must be square to compute its eigenvalues.");
                        }
                    }
                    
                    virtual ~QR() = default;
                    
                    const Matrix<ValueType> &q() {
                        if (_q.is_null()) {
                            ComputeQR(original);
                        }
                        return _q;
                    }
                    
                    const Matrix<ValueType> &r() {
                        if (_r.is_null()) {
                            ComputeQR(original);
                        }
                        return _r;
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
                        if (eigen_values.is_empty()) {
                            auto matrix(original);
                            const size_t last_row = rows - 1;
                            ValueType square_sum, element;
                            
                            for (size_t iteration = 0; iteration < max_iterations; ++iteration) {
                                ComputeQR(matrix);
                                matrix = _r * _q;
                                
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
                            
                            eigen_values = matrix.get_diagonal();
                        }
                        
                        return eigen_values;
                    }
                    
                    const containers::Vector<ValueType> &EigenVector(const ValueType &eigen_value) {
                        
                        auto it_eigen_vector = eigen_vectors.find(eigen_value);
                        if (it_eigen_vector != eigen_vectors.end() && !it_eigen_vector->second.is_empty()) {
                            return it_eigen_vector->second;
                        }
                        
                        Matrix<ValueType> inverse_matrix(rows, rows, 0);
                        inverse_matrix.set_diagonal(eigen_value * (accuracy + 1.0));
                        
                        inverse_matrix = (original - inverse_matrix).pow(-1);
                        
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
                        
                        eigen_vectors.emplace(eigen_value, eigenVector.get_column_as_vector(0) / eigenVector[0][rows - 1]);
                        
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
                    
                    Matrix<ValueType> _q, _r;
                    containers::Vector<ValueType> eigen_values;
                    std::map<ValueType, containers::Vector<ValueType>> eigen_vectors;
                    
                    void ComputeQR(const Matrix<ValueType> &matrix) {
                        
                        const auto I = Matrix<ValueType>::identity(rows);
                        
                        //  First column
                        auto c = matrix.get_column_as_vector(0);
                        auto vt = I.get_row_as_vector(0);
                        vt *= signum(c[0]) * c.norm();
                        vt += c;
                        
                        _q = containers::transpose(vt) * vt;
                        _q *= -2.0 / vt.square_norm();
                        _q += I;
                        
                        _r = _q * matrix;
                        for (auto it_r = _r.begin() + rows; it_r != _r.end(); it_r += rows) {
                            *it_r = 0;
                        }
                        
                        //  Other columns
                        Matrix<ValueType> h, h__;
                        const auto last_row = rows - 1;
                        ValueType *it_r_column, *it_end_r_column;
                        
                        for (size_t row = 1; row < last_row; ++row) {
                            c = _r.get_column_as_vector(row, row);
                            
                            if (c.is_null()) {
                                h = I;
                            } else {
                                vt = I.get_row_as_vector(row, row);
                                vt *= signum(c[0]) * c.norm();
                                vt += c;
                                
                                h = containers::transpose(vt) * vt;
                                h *= -2.0 / vt.square_norm();
                                h += Matrix<ValueType>::identity(rows - row);
                                
                                h__ = I;
                                h__.set_matrix(row, row, h);
                                h = std::move(h__);
                            }
                            
                            _q = h * _q;
                            _r = h * _r;
                            
                            it_r_column = _r[row];
                            it_end_r_column = it_r_column + row;
                            for ( ; it_r_column != it_end_r_column; ++it_r_column) {
                                *it_r_column = 0;
                            }
                        }
                        
                        _q = std::move(_q.transpose());
                        
                        for (size_t row = 1; row < rows; ++row) {
                            it_r_column = _r[row];
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
