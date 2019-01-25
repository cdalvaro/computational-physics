//
//  matrix.hpp
//  Vibration Normal Modes 2D
//
//  Created by Carlos David on 04/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <type_traits>

#include "../algorithms/find.hpp"
#include "../algorithms/factorization/lu.hpp"
#include "vector.hpp"


namespace cda {
    namespace math {
        namespace containers {
        
            //  --- MATRIX CLASS ---
            template <typename T>
            class Matrix {
            private:
                
                typedef typename std::conditional<std::is_floating_point<T>::value, T, float>::type lu_value_type;
                
                size_t n, m, mat_size;
                T *a, *it_end;
                
                void alloc_memory(const size_t &size) {
                    if (size == 0) {
                        std::free(a);
                        a = it_end = nullptr;
                    } else {
                        if (void *mem = std::realloc(a, size * sizeof(T))) {
                            a = static_cast<T *>(mem);
                            it_end = a + size;
                        } else {
                            throw std::bad_alloc();
                        }
                    }
                }
                
            public:
                
                typedef T value_type;
                
                Matrix(const size_t &rows = 0, const size_t &columns = 0) :
                n(rows), m(columns), mat_size(rows * columns),
                a(nullptr), it_end(nullptr) {
                    alloc_memory(mat_size);
                }
                
                Matrix(const size_t &rows, const size_t &columns, const value_type &value) :
                n(rows), m(columns), mat_size(rows * columns),
                a(nullptr), it_end(nullptr) {
                    alloc_memory(mat_size);
                    std::fill(this->begin(), this->end(), value);
                }
                
                Matrix(const Matrix<value_type> &matrix) :
                n(matrix.n), m(matrix.m), mat_size(matrix.mat_size),
                a(nullptr), it_end(nullptr) {
                    alloc_memory(mat_size);
                    std::copy(matrix.begin(), matrix.end(), this->begin());
                }
                
                template<typename ValueType2>
                Matrix(const Matrix<ValueType2> &matrix) :
                n(matrix.rows()), m(matrix.columns()), mat_size(matrix.size()),
                a(nullptr), it_end(nullptr) {
                    alloc_memory(mat_size);
                    std::copy(matrix.begin(), matrix.end(), this->begin());
                }
                
                Matrix(Matrix<value_type> &&matrix) :
                n(matrix.n), m(matrix.m), mat_size(matrix.mat_size),
                a(matrix.a), it_end(matrix.it_end) {
                    matrix.n = matrix.m = matrix.mat_size = 0;
                    matrix.a = matrix.it_end = nullptr;
                }
                
                template <size_t size>
                Matrix(const size_t &rows, const size_t &columns, const value_type (&values)[size]):
                n(rows), m(columns), mat_size(rows * columns),
                a(nullptr), it_end(nullptr) {
                    if (this->mat_size != size) {
                        throw std::logic_error("Sizes do not match");
                    }
                    alloc_memory(size);
                    std::copy(values, values + size, this->begin());
                }
                
                template <size_t rows, size_t columns>
                Matrix(const value_type (&values)[rows][columns]):
                n(rows), m(columns), mat_size(rows * columns),
                a(nullptr), it_end(nullptr) {
                    alloc_memory(mat_size);
                    for (size_t row = 0; row < rows; ++row) {
                        std::copy(values[row], values[row] + columns, this->operator[](row));
                    }
                }
                
                ~Matrix() {
                    if (a) {
                        std::free(a);
                        a = it_end = nullptr;
                        n = m = mat_size = 0;
                    }
                }
                
                value_type *begin() const {
                    return a;
                }
                
                value_type *end() const {
                    return it_end;
                }
                
                void resize(const size_t &rows, const size_t &columns, const bool &fill = false) {
                    if (rows == n && columns == m) {
                        return;
                    }
                    
                    if (columns == m) {
                        const auto new_size = rows * columns;
                        alloc_memory(new_size);
                        
                        if (fill && rows > n) {
                            std::fill_n(a + mat_size, new_size - mat_size, static_cast<value_type>(0));
                        }
                        
                        n = rows;
                        mat_size = new_size;
                    } else {
                        
                        Matrix<value_type> tmp(rows, columns);
                        if (fill) {
                            tmp.zero();
                        }
                        
                        if (mat_size != 0) {
                            const auto &rows_to_copy = rows < n ? rows : n;
                            const auto &columns_to_copy = columns < m ? columns : m;
                            for (size_t row = 0; row < rows_to_copy; ++row) {
                                for (size_t column = 0; column < columns_to_copy; ++column) {
                                    tmp[row][column] = this->operator[](row)[column];
                                }
                            }
                        }
                        
                        this->operator=(std::move(tmp));
                    }
                    
                }
                
                void dimensions(const size_t &rows, const size_t &columns) {
                    if (rows * columns != this->mat_size) {
                        throw std::out_of_range("This method does not resize the matrix, just change the dimensions");
                    }
                    
                    this->n = rows;
                    this->m = columns;
                }
                
                Matrix<value_type> &operator=(const Matrix<value_type> &matrix) {
                    if (this != &matrix) {
                        resize(matrix.n, matrix.m);
                        std::copy(matrix.begin(), matrix.end(), this->begin());
                    }
                    
                    return *this;
                }
                
                Matrix<value_type> &operator=(Matrix<value_type> &&matrix) {
                    if (this != &matrix) {
                        std::free(a);
                        a = matrix.a;
                        it_end = matrix.it_end;
                        n = matrix.n;
                        m = matrix.m;
                        mat_size = matrix.mat_size;
                        
                        matrix.a = matrix.it_end = nullptr;
                        matrix.n = matrix.m = matrix.mat_size = 0;
                    }
                    
                    return *this;
                }
                
                bool operator==(const Matrix<value_type> &matrix) const {
                    if (this->n != matrix.n || this->m != matrix.m) {
                        return false;
                    }
                    
                    auto it_matrix = matrix.begin();
                    for (auto it = this->begin(); it != this->end(); ++it, ++it_matrix) {
                        if (*it != *it_matrix) {
                            return false;
                        }
                    }
                    
                    return true;
                }
                
                bool operator!=(const Matrix<value_type> &matrix) const {
                    return !this->operator==(matrix);
                }
                
                Matrix<value_type> get_row(const size_t &row) const {
                    if (row >= n) {
                        throw std::out_of_range("Index out of bounds");
                    }
                    
                    Matrix<value_type> tmp(1, m);
                    
                    const auto &it_row = this->operator[](row);
                    std::copy(it_row, it_row + m, tmp.begin());
                    
                    return tmp;
                }
                
                Matrix<value_type> get_column(const size_t &column) const {
                    if (column >= m) {
                        throw std::out_of_range("Index out of bounds");
                    }
                    
                    Matrix<value_type> tmp(n, 1);
                    auto it_tmp = tmp.begin();
                    
                    for (size_t i = 0; i < n; ++i, ++it_tmp) {
                        *it_tmp = this->operator[](i)[column];
                    }
                    
                    return tmp;
                }
                
                Matrix<value_type> get_matrix(const size_t &row, const size_t &column,
                                              const size_t &number_of_rows, const size_t &number_of_columns) const {
                    if (row + number_of_rows > n || column + number_of_columns > m) {
                        throw std::out_of_range("Index out of bounds");
                    }
                    
                    Matrix<value_type> tmp(number_of_rows, number_of_columns);
                    
                    auto it_this = begin() + row * m + column;
                    for (auto it_tmp = tmp.begin(); it_tmp != tmp.end(); it_tmp += number_of_columns, it_this += m) {
                        std::copy(it_this, it_this + number_of_columns, it_tmp);
                    }
                    
                    return tmp;
                }
                
                Matrix<value_type> get_matrix(const size_t &row, const size_t &column) const {
                    return this->get_matrix(row, column, n - row, m - column);
                }
                
                Vector<value_type> get_diagonal() const {
                    if (!this->IsSquare()) {
                        throw std::logic_error("Matrix must be an square matrix");
                    }
                    
                    Vector<value_type> diagonal(n);
                    
                    auto it_this = begin();
                    for (auto it_diagonal = diagonal.begin(); it_diagonal != diagonal.end(); ++it_diagonal, it_this += m + 1) {
                        *it_diagonal = *it_this;
                    }
                    
                    return diagonal;
                }
                
                Vector<value_type> get_row_as_vector(const size_t &row, const size_t &from_column = 0) const {
                    if (row >= this->n || from_column >= this->m) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    auto it_begin_row = this->operator[](row) + from_column;
                    auto it_end_row = it_begin_row + this->m - from_column;
                    
                    return Vector<value_type>(it_begin_row, it_end_row);
                }
                
                Vector<value_type> get_column_as_vector(const size_t &column, const size_t &from_row = 0) const {
                    if (column >= this->m || from_row >= this->n) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    Vector<value_type> vector(this->n - from_row);
                    
                    auto it_this = this->begin() + this->m * from_row + column;
                    for (auto it_vector = vector.begin(); it_vector != vector.end(); ++it_vector, it_this += m) {
                        *it_vector = *it_this;
                    }
                    
                    return vector;
                }
                
                //  Sets
                void set_column(const size_t &column, const Matrix<value_type> &matrix) {
                    if (column >= this->m) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    if (matrix.m > 1) {
                        throw std::logic_error("Source matrix must have only one column");
                    }
                    
                    if (matrix.n != this->n) {
                        throw std::logic_error("Matrices must have the same number of rows");
                    }
                    
                    auto it_this = this->begin() + column;
                    for (auto it_matrix = matrix.begin(); it_matrix != matrix.end(); ++it_matrix, it_this += m) {
                        *it_this = *it_matrix;
                    }
                }
                
                void set_column(const size_t &column, const Vector<value_type> &vector) {
                    if (column >= this->m) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    if (vector.size() != this->n) {
                        throw std::logic_error("The vector and the column of the matrix must have the same number of elements");
                    }
                    
                    auto it_this = this->begin() + column;
                    for (auto it_vector = vector.begin(); it_vector != vector.end(); ++it_vector, it_this += m) {
                        *it_this = *it_vector;
                    }
                }
                
                void set_row(const size_t &row, const Matrix<value_type> &matrix) {
                    if (row >= this->n) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    if (matrix.n > 1) {
                        throw std::logic_error("Source matrix must have only one row");
                    }
                    
                    if (matrix.m != this->m) {
                        throw std::logic_error("Matrices must have the same number of columns");
                    }
                    
                    auto it_this = this->operator[](row);
                    std::copy(matrix.begin(), matrix.end(), it_this);
                }
                
                void set_row(const size_t &row, const Vector<value_type> &vector) {
                    if (row >= this->n) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    if (vector.size() != this->m) {
                        throw std::logic_error("The vector and the row of the matrix must have the same number of elements");
                    }
                    
                    auto it_this = this->operator[](row);
                    std::copy(vector.begin(), vector.end(), it_this);
                }
                
                void set_diagonal(const Vector<value_type> &diagonal) {
                    if (!IsSquare()) {
                        throw std::logic_error("The matrix must be square to establish its diagonal");
                    }
                    
                    if (diagonal.size() != this->n) {
                        throw std::logic_error("Diagonal size does not match with the number of rows of the matrix");
                    }
                    
                    auto it_matrix = this->begin();
                    const size_t step = m + 1;
                    for (auto it_diagonal = diagonal.begin(); it_diagonal != diagonal.end(); ++it_diagonal, it_matrix += step) {
                        *it_matrix = *it_diagonal;
                    }
                }
                
                void set_diagonal(const value_type &diagonal) {
                    if (!IsSquare()) {
                        throw std::logic_error("The matrix must be square to establish its diagonal");
                    }
                    
                    auto it_end = end() + m;
                    const size_t step = m + 1;
                    for (auto it = begin(); it != it_end; it += step) {
                        *it = diagonal;
                    }
                }
                
                void set_matrix(const size_t &row, const size_t &column, const Matrix<value_type> &matrix) {
                    if (row >= this->n || column >= this->m) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    const auto number_of_rows = std::min(this->n - row, matrix.rows());
                    const auto number_of_columns = std::min(this->m - column, matrix.columns());
                    
                    auto it_matrix = matrix.begin();
                    const auto it_this_end = this->operator[](row + number_of_rows);
                    for (auto it_this = this->operator[](row); it_this != it_this_end; it_this += this->m, it_matrix += matrix.columns()) {
                        std::copy_n(it_matrix, number_of_columns, it_this + column);
                    }
                }
                
                //  --- FUNCTIONS ---
                //  Returns an integer / a double / a float
                std::pair<size_t, size_t> dimensions() const {
                    return std::pair<size_t, size_t>{ n, m };
                }
                
                size_t rows() const {
                    return n;
                }
                
                size_t columns() const {
                    return m;
                }
                
                size_t size() const {
                    return mat_size;
                }
                
                Matrix<value_type> sum_rows() const {
                    Matrix<value_type> sum_rows(n, 1, 0);
                    
                    auto it_this = this->begin();
                    for (auto it_sum = sum_rows.begin(); it_sum != sum_rows.end(); ++it_sum) {
                        auto it_this_row_end = it_this + this->m;
                        for (; it_this != it_this_row_end; ++it_this) {
                            *it_sum += *it_this;
                        }
                    }
                    
                    return sum_rows;
                }
                
                Vector<value_type> sum_rows_as_vector() const {
                    return sum_rows().get_column_as_vector(0);
                }
                
                value_type sum_row(const size_t &row) const {
                    if (row >= this->n) {
                        throw std::out_of_range("Index out of range");
                    }
                    
                    value_type sum = 0;
                    
                    const auto &it_end = this->operator[](row) + this->m;
                    for (auto it_row = this->operator[](row); it_row != it_end; ++it_row) {
                        sum += *it_row;
                    }
                    
                    return sum;
                }
                
                Matrix<value_type> sum_columns() const {
                    Matrix<value_type> sum_columns(1, m, 0);
                    
                    for (auto it_this = this->begin(); it_this != this->end();) {
                        for (auto it_sum = sum_columns.begin(); it_sum != sum_columns.end(); ++it_sum, ++it_this) {
                            *it_sum += *it_this;
                        }
                    }
                    
                    return sum_columns;
                }
                
                Vector<value_type> sum_columns_as_vector() const {
                    return sum_columns().get_row_as_vector(0);
                }
                
                value_type SumColumn(const size_t &column) const {
                    if (column >= this->m) {
                        throw std::out_of_range("Index out of range");
                    }
                    
                    value_type sum = 0;
                    
                    const auto &it_column_end = this->end() + column;
                    for (auto it_column = this->begin() + column; it_column != it_column_end; it_column += this->m) {
                        sum += *it_column;
                    }
                    
                    return sum;
                }
                
                value_type max_element() const {
                    return algorithms::find::max_element(begin(), end());
                }
                
                value_type abs_max_element() const {
                    return algorithms::find::abs_max_element(begin(), end());
                }
                
                value_type abs_max_element_with_sign() const {
                    return algorithms::find::abs_max_element_with_sign(begin(), end());
                }
                
                value_type min_element() const {
                    return algorithms::find::min_element(begin(), end());
                }
                
                value_type abs_min_element() const {
                    return algorithms::find::abs_min_element(begin(), end());
                }
                
                value_type abs_min_element_with_sign() const {
                    return algorithms::find::abs_min_element_with_sign(begin(), end());
                }
                
                const value_type &at(const size_t &row, const size_t &column) const {
                    if (row >= n || column >= m) {
                        throw std::out_of_range("Indexes out of range");
                    }
                    return this->operator[](row)[column];
                }
                
                const value_type *operator[](const size_t &row) const {
                    return &a[row * m];
                }
                
                value_type &at(const size_t &row, const size_t &column) {
                    if (row >= n || column >= m) {
                        throw std::out_of_range("Indexes out of range");
                    }
                    return this->operator[](row)[column];
                }
                
                value_type *operator[](const size_t &row) {
                    return &a[row * m];
                }
                
                Matrix<value_type> Transpose() const {
                    Matrix<value_type> new_matrix(this->m, this->n);
                    
                    for (size_t row = 0; row < this->n; ++row) {
                        for (size_t column = 0; column < this->m; ++column) {
                            new_matrix[column][row] = this->operator[](row)[column];
                        }
                    }
                    
                    return new_matrix;
                }
                
                void clear() {
                    alloc_memory(0);
                    n = m = mat_size = 0;
                }
                
                bool is_empty() const {
                    return mat_size == 0;
                }
                
                bool is_null() const {
                    for (auto it = begin(); it != end(); ++it) {
                        if (*it != 0) {
                            return false;
                        }
                    }
                    return true;
                }
                
                bool IsSquare() const {
                    return n == m;
                }
                
                bool has_duplicate(const value_type &accuracy) const {
                    value_type distance;
                    for (auto it1 = begin(); it1 != end(); ++it1) {
                        for (auto it2 = begin(); it2 != end(); ++it2) {
                            if (it1 != it2) {
                                distance = *it1 - *it2;
                                if (std::sqrt(distance * distance) < accuracy) {
                                    return true;
                                }
                            }
                        }
                    }
                    return false;
                }
                
                bool has_duplicate() const {
                    for (auto it1 = begin(); it1 != end(); ++it1) {
                        for (auto it2 = begin(); it2 != end(); ++it2) {
                            if (it1 != it2 && *it1 == *it2) {
                                return true;
                            }
                        }
                    }
                    return false;
                }
                
                value_type * find(const value_type &value) const {
                    return cda::math::algorithms::find::Element(begin(), end(), value);
                }
                
                //  Void functions
                void fill(const value_type &value) {
                    std::fill(begin(), end(), value);
                }
                
                void zero() {
                    fill(0);
                }
                
                void ones() {
                    fill(1);
                }
                
                void Identity() {
                    if (!IsSquare()) {
                        throw std::logic_error("Unable to create an identity matrix in a non squere matrix.");
                    }
                    
                    zero();
                    auto it_end = end() + m;
                    const size_t step = m + 1;
                    for (auto it = begin(); it != it_end; it += step) {
                        *it = 1;
                    }
                }
                
                Matrix<value_type> operator+(const Matrix<value_type> &matrix) const {
                    if (this->n != matrix.n || this->m != matrix.m) {
                        throw std::logic_error("Matrices must be of the same dimensions.");
                    }
                    
                    Matrix<value_type> new_matrix(n, m);
                    auto it_new = new_matrix.begin();
                    auto it_matrix = matrix.begin();
                    
                    for (auto it_this = this->begin(); it_this != this->end(); ++it_this, ++it_matrix, ++it_new) {
                        *it_new = *it_this + *it_matrix;
                    }
                    
                    return new_matrix;
                }
                
                Matrix<value_type>& operator+=(const Matrix<value_type> &matrix) {
                    if (this->n != matrix.n || this->m != matrix.m) {
                        throw std::logic_error("Matrices must be of the same dimensions.");
                    }
                    
                    auto it_this = this->begin();
                    for (auto it_matrix = matrix.begin(); it_matrix != matrix.end(); ++it_matrix) {
                        *it_this++ += *it_matrix;
                    }
                    
                    return *this;
                }
                
                Matrix<value_type> operator-(const Matrix<value_type> &matrix) const {
                    if (this->n != matrix.n || this->m != matrix.m) {
                        throw std::logic_error("Matrices must be of the same dimensions.");
                    }
                    
                    Matrix<value_type> new_matrix(n, m);
                    auto it_new = new_matrix.begin();
                    auto it_matrix = matrix.begin();
                    
                    for (auto it_this = this->begin(); it_this != this->end(); ++it_this, ++it_matrix, ++it_new) {
                        *it_new = *it_this - *it_matrix;
                    }
                    
                    return new_matrix;
                }
                
                Matrix<value_type>& operator-=(const Matrix<value_type> &matrix) {
                    if (this->n != matrix.n || this->m != matrix.m) {
                        throw std::logic_error("Matrices must be of the same dimensions.");
                    }
                    
                    auto it_this = this->begin();
                    for (auto it_matrix = matrix.begin(); it_matrix != matrix.end(); ++it_matrix) {
                        *it_this++ -= *it_matrix;
                    }
                    
                    return *this;
                }
                
                Matrix<value_type> operator*(const value_type &value) const {
                    Matrix<value_type> new_matrix(n, m);
                    auto it_new = new_matrix.begin();
                    
                    for (auto it_this = this->begin(); it_this != this->end(); ++it_this) {
                        *it_new++ = *it_this * value;
                    }
                    
                    return new_matrix;
                }
                
                Matrix<value_type> operator*(const Matrix<value_type> &matrix) const {
                    if (this->m != matrix.n) {
                        throw std::logic_error("Matrices dimensions are not compatible.");
                    }
                    
                    Matrix<value_type> new_matrix(this->n, matrix.m);
                    auto it_new_matrix = new_matrix.begin();
                    
                    const value_type *it_this_row, *it_this_row__, *it_end_this_row;
                    
                    value_type sum;
                    size_t matrix_column;
                    value_type *it_matrix_column;
                    
                    for (it_this_row = this->begin(); it_this_row < this->end(); it_this_row = it_end_this_row) {
                        it_end_this_row = it_this_row + this->m;
                        
                        for (matrix_column = 0; matrix_column < matrix.m; ++matrix_column) {
                            sum = 0;
                            for (it_this_row__ = it_this_row, it_matrix_column = matrix.begin() + matrix_column;
                                 it_this_row__ != it_end_this_row;
                                 ++it_this_row__, it_matrix_column += matrix.m) {
                                sum += *it_this_row__ * *it_matrix_column;
                            }
                            *it_new_matrix++ = sum;
                        }
                    }
                    
                    return new_matrix;
                }
                
                Matrix<value_type> &operator*=(const Matrix<value_type> &matrix) {
                    *this = *this * matrix;
                    return *this;
                }
                
                Matrix<value_type> operator*(const Vector<value_type> &vector) {
                    if (m != 1) {
                        throw std::logic_error("Matrix and Vector are not compatible.");
                    }
                    
                    Matrix<value_type> new_matrix(n, vector.size());
                    auto it_new = new_matrix.begin();
                    
                    for (auto it_this = this->begin(); it_this != this->end(); ++it_this) {
                        for (auto it_vector = vector.begin(); it_vector != vector.end(); ++it_vector) {
                            *it_new = *it_this * *it_vector;
                            ++it_new;
                        }
                    }
                    
                    return new_matrix;
                }
                
                Matrix<value_type>& operator*=(const value_type &value) {
                    for (auto it = begin(); it != end(); ++it) {
                        *it *= value;
                    }
                    return *this;
                }
                
                Matrix<value_type> operator/(const value_type &value) const {
                    Matrix<value_type> new_matrix(this->n, this->m);
                    auto it_new = new_matrix.begin();
                    
                    for (auto it_this = this->begin(); it_this != this->end(); ++it_this, ++it_new) {
                        *it_new = *it_this / value;
                    }
                    
                    return new_matrix;
                }
                
                Matrix<value_type>& operator/=(const value_type &value) {
                    for (auto it = begin(); it != end(); ++it) {
                        *it /= value;
                    }
                    return *this;
                }
                
                Matrix<value_type> operator-() const {
                    Matrix<value_type> new_matrix(this->n, this->m);
                    std::transform(this->begin(), this->end(), new_matrix.begin(),
                                   [](const value_type &value) {
                                       return -value;
                                   });
                    
                    return new_matrix;
                }
                
                value_type Determinant() const {
                    return algorithms::factorization::LU<Matrix, lu_value_type>::Determinant(*this);
                }
                
                Matrix<value_type> Pow(const ssize_t &power) const {
                    if (!this->IsSquare()) {
                        throw std::logic_error("Matrix must be an square matrix");
                    }
                    
                    Matrix<value_type> new_matrix;
                    
                    switch (power) {
                        case -1:
                            new_matrix = algorithms::factorization::LU<Matrix, lu_value_type>::InverseMatrix(*this);
                            break;
                            
                        case 0:
                            if (Determinant() == 0) {
                                throw std::logic_error("Matrix is singular");
                            }
                            new_matrix = Identity(this->n);
                            break;
                            
                        default:
                            new_matrix = *this;
                            for (size_t i = 2; i <= power; ++i) {
                                new_matrix *= *this;
                            }
                            break;
                    }
                    
                    return new_matrix;
                }
                
                //  --- STATIC METHODS ---
                static Matrix<value_type> zero(const size_t &rows, const size_t &columns) {
                    return Matrix<value_type>(rows, columns, 0);
                }
                
                static Matrix<value_type> ones(const size_t &rows, const size_t &columns) {
                    return Matrix<value_type>(rows, columns, 1);
                }
                
                static Matrix<value_type> Identity(const size_t &rows) {
                    Matrix<value_type> identity(rows, rows);
                    identity.Identity();
                    return identity;
                }
                
            };
            
            //  MARK: - Extra funcions
            template <typename ValueType>
            Matrix<ValueType> Transpose(const Vector<ValueType> &vector) {
                Matrix<ValueType> matrix(vector.size(), 1);
                std::copy(vector.begin(), vector.end(), matrix.begin());
                return matrix;
            }
        
        } /* namespace math */
    } /* namespace containers */
} /* namespace cda */

//  MARK: - Extra operators
template <typename ValueType>
cda::math::containers::Matrix<ValueType> operator*(const ValueType &value,
                                           const cda::math::containers::Matrix<ValueType> &matrix) {
    
    cda::math::containers::Matrix<ValueType> tmp(matrix.rows(), matrix.columns());
    auto it_tmp = tmp.begin();
    
    for (auto it_matrix = matrix.begin(); it_matrix != matrix.end(); ++it_matrix, ++it_tmp) {
        *it_tmp = *it_matrix * value;
    }
    
    return tmp;
}

template <typename ValueType>
cda::math::containers::Vector<ValueType> operator*(const cda::math::containers::Vector<ValueType> &vector,
                                           const cda::math::containers::Matrix<ValueType> &matrix) {
    
    const auto &rows = matrix.rows();
    if (rows != vector.size()) {
        throw std::logic_error("The vector and the matrix are incompatible");
    }
    
    const auto &columns = matrix.columns();
    cda::math::containers::Vector<ValueType> new_vector(columns, 0);
    
    auto it_new = new_vector.begin();
    for (size_t column = 0; column < columns; ++column, ++it_new) {
        auto it_vector = vector.begin();
        for (size_t row = 0; row < rows; ++row, ++it_vector) {
            *it_new += *it_vector * matrix[row][column];
        }
    }
    
    return new_vector;
}

template <typename ValueType>
cda::math::containers::Matrix<ValueType> operator&&(const cda::math::containers::Matrix<ValueType> &left_matrix,
                                            const cda::math::containers::Matrix<ValueType> &right_matrix) {
    
    if (left_matrix.rows() != right_matrix.rows()) {
        throw std::logic_error("Both matrices must have the same number of rows");
    }
    
    cda::math::containers::Matrix<ValueType> new_matrix(left_matrix.rows(), left_matrix.columns() + right_matrix.columns());
    new_matrix.set_matrix(0, 0, left_matrix);
    new_matrix.set_matrix(0, left_matrix.columns(), right_matrix);
    
    return new_matrix;
}

template <typename ValueType>
void operator||(cda::math::containers::Matrix<ValueType> &left_matrix,
                cda::math::containers::Matrix<ValueType> &right_matrix) {
    
    const size_t left_matrix_columns = left_matrix.columns() / 2;
    
    right_matrix = left_matrix.get_matrix(0, left_matrix_columns, left_matrix.rows(), left_matrix.columns() - left_matrix_columns);
    left_matrix.resize(left_matrix.rows(), left_matrix_columns);
}

template <typename ValueType>
void operator>>(std::istream &input,
                cda::math::containers::Matrix<ValueType> &matrix) {
    
    if (!input) {
        throw std::logic_error("Input is not avaiable");
    }
    
    size_t initial_size = 100;
    matrix.resize(initial_size, 1);
    
    size_t rows = 0, columns = 0;
    
    bool is_first_row = true;
    char comma;
    std::string line, cell;
    size_t element = 0;
    while (std::getline(input, line)) {
        
        std::stringstream line_stream(line);
        while (line_stream >> matrix[element][0]) {
            
            ++element;
            
            if (element >= matrix.rows()) {
                initial_size *= 2;
                matrix.resize(matrix.rows() + initial_size, 1);
            }
            
            if (is_first_row) {
                ++columns;
            }
            
            line_stream >> comma;
        }
        
        if (is_first_row) {
            is_first_row = false;
        }
        
        ++rows;
    }
    
    if (rows * columns != element) {
        throw std::out_of_range("The size of the matrix could not be determined properly.");
    }
    
    matrix.resize(element, 1);
    matrix.dimensions(rows, columns);
}

template <typename ValueType>
std::ostream& operator<<(std::ostream &output,
                         const cda::math::containers::Matrix<ValueType> &matrix) {
    
    const size_t rows = matrix.rows();
    const size_t columns = matrix.columns();
    
    if (output.rdbuf() == std::cout.rdbuf()) {
        
        const size_t custom_width = 12;
        const size_t custom_precision = 5;
        
        output.width();
        output << std::fixed;
        output.fill(' ');
        output.precision(custom_precision);
    
        if (rows == 1 && columns == 1) {
            output << matrix[0][0] << std::endl;
        } else if (rows == 1) {
            output << "[";
            output.width(custom_width);
            output << matrix[0][0];
            for (size_t column = 1; column < columns; ++column) {
                output << " ";
                output.width(custom_width);
                output << matrix[0][column];
            }
            output << "]" << std::endl;
        } else {
            
            const size_t last_row = rows - 1;
            
            output << "⎡";
            output.width(custom_width);
            output << matrix[0][0];
            for (size_t column = 1; column < columns; ++column) {
                output << std::right << " ";
                output.width(custom_width);
                output << matrix[0][column];
            }
            output.width();
            output << " ⎤" << std::endl;
            for (size_t row = 1; row < last_row; ++row) {
                output << "⎢";
                output.width(custom_width);
                output << matrix[row][0];
                for (size_t column = 1; column < columns; ++column) {
                    output << std::right << " ";
                    output.width(custom_width);
                    output << matrix[row][column];
                }
                output << " ⎥" << std::endl;
            }
            output << "⎣";
            output.width(custom_width);
            output << matrix[last_row][0];
            for (size_t column = 1; column < columns; ++column) {
                output << std::right << " ";
                output.width(custom_width);
                output << matrix[last_row][column];
            }
            output << " ⎦" << std::endl;
        }
        output.precision();
    
    } else {
        const char separator = ';';
        for (size_t row = 0; row < rows; ++row) {
            auto it_row = matrix[row];
            auto it_end_row = it_row + columns;
            
            output << *it_row;
            for (it_row = std::next(it_row); it_row != it_end_row; ++it_row) {
                output << separator << *it_row;
            }
            output << std::endl;
        }
    }

    return output;
}
