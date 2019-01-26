//
//  linear.hpp
//  Computational Physics
//
//  Created by Carlos David on 17/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#pragma once

#include "../../containers.hpp"
#include "../../algorithms/factorization/lu.hpp"


#define CDA_LINEAR_DEFAULT_ACCURACY 1E-06

namespace cda {
    namespace math {
        namespace equations {
            namespace systems {
                namespace linear {
                    
                    template <typename T>
                    containers::Vector<T> solve_lu(const containers::Matrix<T> &system,
                                                   const containers::Vector<T> &b_terms) {
                        algorithms::factorization::LU<containers::Matrix, T> lu(system);
                        return lu.solve_linear_system(b_terms);
                    }
                    
                    template <typename T>
                    containers::Vector<T> solve_3diagonal(const containers::Matrix<T> &system,
                                                          const containers::Vector<T> &b_terms) {
                        
                        if (!system.is_square()) {
                            throw std::logic_error("The system is matrix is not square");
                        }
                        
                        auto rows = system.rows();
                        
                        if (rows != b_terms.size()) {
                            throw std::logic_error("The number of rows of the system matrix does not match the number of elements in the b terms vector.");
                        }
                        
                        containers::Vector<T> tmp(rows, 0), x(rows, 0), al(rows, 0), be(rows, 0);
                        
                        be[0] = system[0][1];
                        for (size_t row = 1; row < rows; ++row) {
                            al[row] = system[row][0] / be[row - 1];
                            be[row] = system[row][1] - al[row] * system[row - 1][2];
                        }
                        
                        tmp[0] = b_terms[0];
                        for (size_t row = 1; row < rows; ++row) {
                            tmp[row] = b_terms[row] - al[row] * tmp[row - 1];
                        }
                        
                        x[rows - 1] = tmp[rows - 1] / be[rows - 1];
                        for (ssize_t row = rows - 2; row >= 0; --row) {
                            x[row] = (tmp[row] - system[row][2] * x[row + 1]) / be[row];
                        }
                        
                        return x;
                    }
                    
                    template<typename T>
                    containers::Vector<T> solve_gauss_seidel_3diagonal(const containers::Matrix<T> &system,
                                                                       const containers::Vector<T> &b_terms,
                                                                       const double &accuracy = CDA_LINEAR_DEFAULT_ACCURACY) {
                        
                        if (!system.is_square()) {
                            throw std::logic_error("The system is matrix is not square");
                        }
                        
                        const auto rows = system.rows();
                        
                        if (rows != b_terms.size()) {
                            throw std::logic_error("The number of rows of the system matrix does not match the number of elements in the b terms vector.");
                        }
                        
                        containers::Vector<T> x(rows, 0);
                        T error, prev_x, aux_calc;
                        do {
                            error = 0.0;
                            
                            prev_x = x[0];
                            x[0] = (b_terms[0] - system[0][0] * x[1]) / system[0][1];
                            
                            aux_calc = x[0] - prev_x;
                            error += aux_calc * aux_calc;
                            
                            for (size_t row = 1; row < rows - 1; ++row) {
                                prev_x = x[row];
                                x[row] = (b_terms[row] - system[row][0] * x[row - 1] - system[row][2] * x[row + 1]) / system[row][1];
                                
                                aux_calc = x[row] - prev_x;
                                error += aux_calc * aux_calc;
                            }
                            
                            prev_x = x[rows-1];
                            x[rows - 1] = (b_terms[rows - 1] - system[rows - 1][0] * x[rows-2]) / system[rows - 1][1];
                            
                            aux_calc = x[rows - 1] - prev_x;
                            error += aux_calc * aux_calc;
                            
                            error = std::sqrt(error);
                            
                        } while (error > accuracy);
                        
                        return x;
                    }
                    
                }
            }
        }
    }
}
