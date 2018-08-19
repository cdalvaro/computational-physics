//
//  matrix.hpp
//  Vibration Normal Modes 2D
//
//  Created by Carlos David on 04/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#pragma once

#include "vector.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <istream>
#include <sstream>
#include <stdlib.h>
#include <string>


namespace cda {
    namespace math {
        namespace containers {
        
            //  --- MATRIX CLASS ---
            template <typename T>
            class Matrix {
            private:
                size_t n, m, size;
                T *a, *it_end;
                
                void AllocateMemory(const size_t &size) {
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
                
                Matrix(const size_t &rows = 0, const size_t &columns = 0) :
                n(rows), m(columns), size(rows * columns),
                a(nullptr), it_end(nullptr) {
                    AllocateMemory(size);
                }
                
                Matrix(const size_t &rows, const size_t &columns, const T &value) :
                n(rows), m(columns), size(rows * columns),
                a(nullptr), it_end(nullptr) {
                    AllocateMemory(size);
                    std::fill(this->Begin(), this->End(), value);
                }
                
                Matrix(const Matrix<T> &matrix) :
                n(matrix.n), m(matrix.m), size(matrix.size),
                a(nullptr), it_end(nullptr) {
                    AllocateMemory(size);
                    std::copy(matrix.Begin(), matrix.End(), this->Begin());
                }
                
                Matrix(Matrix<T> &&matrix) :
                n(matrix.n), m(matrix.m), size(matrix.size),
                a(matrix.a), it_end(matrix.it_end) {
                    matrix.n = matrix.m = matrix.size = 0;
                    matrix.a = matrix.it_end = nullptr;
                }
                
                template <size_t size>
                Matrix(const size_t &rows, const size_t &columns, const T (&values)[size]):
                n(rows), m(columns), size(rows * columns),
                a(nullptr), it_end(nullptr) {
                    if (this->size != size) {
                        throw std::logic_error("Sizes do not match");
                    }
                    AllocateMemory(size);
                    std::copy(values, values + size, this->Begin());
                }
                
                template <size_t rows, size_t columns>
                Matrix(const T (&values)[rows][columns]):
                n(rows), m(columns), size(rows * columns),
                a(nullptr), it_end(nullptr) {
                    AllocateMemory(size);
                    for (size_t row = 0; row < rows; ++row) {
                        std::copy(values[row], values[row] + columns, this->operator[](row));
                    }
                }
                
                ~Matrix() {
                    if (a) {
                        std::free(a);
                        a = it_end = nullptr;
                        n = m = size = 0;
                    }
                }
                
                T *Begin() const {
                    return a;
                }
                
                T *End() const {
                    return it_end;
                }
                
                void Resize(const size_t &rows, const size_t &columns, const bool &fill = false) {
                    if (rows == n && columns == m) {
                        return;
                    }
                    
                    if (columns == m) {
                        const auto new_size = rows * columns;
                        AllocateMemory(new_size);
                        
                        if (fill && rows > n) {
                            std::fill_n(a + size, new_size - size, static_cast<T>(0));
                        }
                        
                        n = rows;
                        size = new_size;
                    } else {
                        
                        Matrix<T> tmp(rows, columns);
                        if (fill) {
                            tmp.Zero();
                        }
                        
                        if (size != 0) {
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
                
                void ChangeDimensions(const size_t &rows, const size_t &columns) {
                    if (rows * columns != this->size) {
                        throw std::out_of_range("This method does not resize the matrix, just change the dimensions");
                    }
                    
                    this->n = rows;
                    this->m = columns;
                }
                
                Matrix<T> &operator=(const Matrix<T> &matrix) {
                    if (this != &matrix) {
                        Resize(matrix.n, matrix.m);
                        std::copy(matrix.Begin(), matrix.End(), this->Begin());
                    }
                    
                    return *this;
                }
                
                Matrix<T> &operator=(Matrix<T> &&matrix) {
                    if (this != &matrix) {
                        std::free(a);
                        a = matrix.a;
                        it_end = matrix.it_end;
                        n = matrix.n;
                        m = matrix.m;
                        size = matrix.size;
                        
                        matrix.a = matrix.it_end = nullptr;
                        matrix.n = matrix.m = matrix.size = 0;
                    }
                    
                    return *this;
                }
                
                bool operator==(const Matrix<T> &matrix) const {
                    if (this->n != matrix.n || this->m != matrix.m) {
                        return false;
                    }
                    
                    auto it_matrix = matrix.Begin();
                    for (auto it = this->Begin(); it != this->End(); ++it, ++it_matrix) {
                        if (*it != *it_matrix) {
                            return false;
                        }
                    }
                    
                    return true;
                }
                
                bool operator!=(const Matrix<T> &matrix) const {
                    return !this->operator==(matrix);
                }
                
                Matrix<T> GetRow(const size_t &row) const {
                    if (row >= n) {
                        throw std::out_of_range("Index out of bounds");
                    }
                    
                    Matrix<T> tmp(1, m);
                    
                    const auto &it_row = this->operator[](row);
                    std::copy(it_row, it_row + m, tmp.Begin());
                    
                    return tmp;
                }
                
                Matrix<T> GetColumn(const size_t &column) const {
                    if (column >= m) {
                        throw std::out_of_range("Index out of bounds");
                    }
                    
                    Matrix<T> tmp(n, 1);
                    auto it_tmp = tmp.Begin();
                    
                    for (size_t i = 0; i < n; ++i, ++it_tmp) {
                        *it_tmp = this->operator[](i)[column];
                    }
                    
                    return tmp;
                }
                
                Matrix<T> GetMatrix(const size_t &row, const size_t &column,
                                    const size_t &number_of_rows, const size_t &number_of_columns) const {
                    if (row + number_of_rows > n || column + number_of_columns > m) {
                        throw std::out_of_range("Index out of bounds");
                    }
                    
                    Matrix<T> tmp(number_of_rows, number_of_columns);
                    
                    auto it_this = Begin() + row * m + column;
                    for (auto it_tmp = tmp.Begin(); it_tmp != tmp.End(); it_tmp += number_of_columns, it_this += m) {
                        std::copy(it_this, it_this + number_of_columns, it_tmp);
                    }
                    
                    return tmp;
                }
                
                Matrix<T> GetMatrix(const size_t &row, const size_t &column) const {
                    return this->GetMatrix(row, column, n - row, m - column);
                }
                
                Vector<T> GetDiagonal() const {
                    if (!this->IsSquared()) {
                        throw std::logic_error("Matrix must be an square matrix");
                    }
                    
                    Vector<T> diagonal(n);
                    
                    auto it_this = Begin();
                    for (auto it_diagonal = diagonal.Begin(); it_diagonal != diagonal.End(); ++it_diagonal, it_this += m + 1) {
                        *it_diagonal = *it_this;
                    }
                    
                    return diagonal;
                }
                
                Vector<T> GetRowAsVector(const size_t &row, const size_t &from_column = 0) const {
                    if (row >= this->n || from_column >= this->m) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    Vector<T> vector(this->m - from_column);
                    vector.Copy(vector.Size(), this->operator[](row) + from_column);
                    
                    return vector;
                }
                
                Vector<T> GetColumnAsVector(const size_t &column, const size_t &from_row = 0) const {
                    if (column >= this->m || from_row >= this->n) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    Vector<T> vector(this->n - from_row);
                    
                    auto it_this = this->Begin() + this->m * from_row + column;
                    for (auto it_vector = vector.Begin(); it_vector != vector.End(); ++it_vector, it_this += m) {
                        *it_vector = *it_this;
                    }
                    
                    return vector;
                }
                
                //  Sets
                void SetColumn(const size_t &column, const Matrix<T> &matrix) {
                    if (column >= this->m) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    if (matrix.m > 1) {
                        throw std::logic_error("Source matrix must have only one column");
                    }
                    
                    if (matrix.n != this->n) {
                        throw std::logic_error("Matrices must have the same number of rows");
                    }
                    
                    auto it_this = this->Begin() + column;
                    for (auto it_matrix = matrix.Begin(); it_matrix != matrix.End(); ++it_matrix, it_this += m) {
                        *it_this = *it_matrix;
                    }
                }
                
                void SetColumn(const size_t &column, const Vector<T> &vector) {
                    if (column >= this->m) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    if (vector.Size() != this->n) {
                        throw std::logic_error("The vector and the column of the matrix must have the same number of elements");
                    }
                    
                    auto it_this = this->Begin() + column;
                    for (auto it_vector = vector.Begin(); it_vector != vector.End(); ++it_vector, it_this += m) {
                        *it_this = *it_vector;
                    }
                }
                
                void SetRow(const size_t &row, const Matrix<T> &matrix) {
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
                    std::copy(matrix.Begin(), matrix.End(), it_this);
                }
                
                void SetRow(const size_t &row, const Vector<T> &vector) {
                    if (row >= this->n) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    if (vector.Size() != this->m) {
                        throw std::logic_error("The vector and the row of the matrix must have the same number of elements");
                    }
                    
                    auto it_this = this->operator[](row);
                    std::copy(vector.Begin(), vector.End(), it_this);
                }
                
                void SetMatrix(const size_t &row, const size_t &column, const Matrix<T> &matrix) {
                    if (row >= this->n || column >= this->m) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    const auto number_of_rows = std::min(this->n - row, matrix.Rows());
                    const auto number_of_columns = std::min(this->m - column, matrix.Columns());
                    
                    auto it_matrix = matrix.Begin();
                    const auto it_this_end = this->operator[](row + number_of_rows);
                    for (auto it_this = this->operator[](row); it_this != it_this_end; it_this += this->m, it_matrix += matrix.Columns()) {
                        std::copy_n(it_matrix, number_of_columns, it_this + column);
                    }
                }
                
                //  --- FUNCTIONS ---
                //  Returns an integer / a double / a float
                std::pair<size_t, size_t> Dimensions() const {
                    return std::pair<size_t, size_t>{ n, m };
                }
                
                size_t Rows() const {
                    return n;
                }
                
                size_t Columns() const {
                    return m;
                }
                
                size_t Size() const {
                    return size;
                }
                
                Matrix<T> SumRows() const {
                    Matrix<T> sum_rows(n, 1, 0);
                    
                    auto it_this = this->Begin();
                    for (auto it_sum = sum_rows.Begin(); it_sum != sum_rows.End(); ++it_sum) {
                        auto it_this_row_end = it_this + this->m;
                        for (; it_this != it_this_row_end; ++it_this) {
                            *it_sum += *it_this;
                        }
                    }
                    
                    return sum_rows;
                }
                
                Vector<T> SumRowsAsVector() const {
                    return SumRows().GetColumnAsVector(0);
                }
                
                T SumRow(const size_t &row) const {
                    if (row >= this->n) {
                        throw std::out_of_range("Index out of range");
                    }
                    
                    T sum = 0;
                    
                    const auto &it_end = this->operator[](row) + this->m;
                    for (auto it_row = this->operator[](row); it_row != it_end; ++it_row) {
                        sum += *it_row;
                    }
                    
                    return sum;
                }
                
                Matrix<T> SumColumns() const {
                    Matrix<T> sum_columns(1, m, 0);
                    
                    for (auto it_this = this->Begin(); it_this != this->End();) {
                        for (auto it_sum = sum_columns.Begin(); it_sum != sum_columns.End(); ++it_sum, ++it_this) {
                            *it_sum += *it_this;
                        }
                    }
                    
                    return sum_columns;
                }
                
                Vector<T> SumColumnsAsVector() const {
                    return SumColumns().GetRowAsVector(0);
                }
                
                T SumColumn(const size_t &column) const {
                    if (column >= this->m) {
                        throw std::out_of_range("Index out of range");
                    }
                    
                    T sum = 0;
                    
                    const auto &it_column_end = this->End() + column;
                    for (auto it_column = this->Begin() + column; it_column != it_column_end; it_column += this->m) {
                        sum += *it_column;
                    }
                    
                    return sum;
                }
                
                T MaximumElement() const {
                    auto it = this->Begin();
                    
                    T max_element = *it;
                    ++it;
                    
                    for (; it != this->End(); ++it) {
                        if (*it > max_element) {
                            max_element = *it;
                        }
                    }
                    
                    return max_element;
                }
                
                T AbsoluteMaximumElement() const {
                    auto it = this->Begin();
                    
                    T max_element = std::abs(*it);
                    ++it;
                    
                    T abs_it;
                    for (; it != this->End(); ++it) {
                        abs_it = std::abs(*it);
                        if (abs_it > max_element) {
                            max_element = abs_it;
                        }
                    }
                    
                    return max_element;
                }
                
                T AbsoluteMaximumElementWithSign() const {
                    auto it = this->Begin();
                    
                    T max_element = *it;
                    ++it;
                    
                    for (; it != this->End(); ++it) {
                        if (std::abs(*it) > std::abs(max_element)) {
                            max_element = *it;
                        }
                    }
                    
                    return max_element;
                }
                
                T MinimumElement() const {
                    auto it = this->Begin();
                    
                    T min_element = *it;
                    ++it;
                    
                    for (; it != this->End(); ++it) {
                        if (*it < min_element) {
                            min_element = *it;
                        }
                    }
                    
                    return min_element;
                }
                
                T AbsoluteMinimumElement() const {
                    auto it = this->Begin();
                    
                    T min_element = std::abs(*it);
                    ++it;
                    
                    T abs_it;
                    for (; it != this->End(); ++it) {
                        abs_it = std::abs(*it);
                        if (abs_it < min_element) {
                            min_element = abs_it;
                        }
                    }
                    
                    return min_element;
                }
                
                T AbsoluteMinimumElementWithSign() const {
                    auto it = this->Begin();
                    
                    T min_element = *it;
                    ++it;
                    
                    for (; it != this->End(); ++it) {
                        if (std::abs(*it) < std::abs(min_element)) {
                            min_element = *it;
                        }
                    }
                    
                    return min_element;
                }
                
                const T &At(const size_t &row, const size_t &column) const {
                    if (row >= n || column >= m) {
                        throw std::out_of_range("Indexes out of range");
                    }
                    return this->operator[](row)[column];
                }
                
                const T *operator[](const size_t &row) const {
                    return &a[row * m];
                }
                
                T &At(const size_t &row, const size_t &column) {
                    if (row >= n || column >= m) {
                        throw std::out_of_range("Indexes out of range");
                    }
                    return this->operator[](row)[column];
                }
                
                T *operator[](const size_t &row) {
                    return &a[row * m];
                }
                
                Matrix<T> Transpose() const {
                    Matrix<T> new_matrix(this->m, this->n);
                    
                    for (size_t row = 0; row < this->n; ++row) {
                        for (size_t column = 0; column < this->m; ++column) {
                            new_matrix[column][row] = this->operator[](row)[column];
                        }
                    }
                    
                    return new_matrix;
                }
                
                void Clear() {
                    AllocateMemory(0);
                    n = m = size = 0;
                }
                
                bool IsEmpty() const {
                    return size == 0;
                }
                
                bool IsNull() const {
                    for (auto it = Begin(); it != End(); ++it) {
                        if (*it != 0) {
                            return false;
                        }
                    }
                    return true;
                }
                
                bool IsSquared() const {
                    return n == m;
                }
                
                bool HasDuplicate(const T &accuracy) const {
                    T distance;
                    for (auto it1 = Begin(); it1 != End(); ++it1) {
                        for (auto it2 = Begin(); it2 != End(); ++it2) {
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
                
                bool HasDuplicate() const {
                    for (auto it1 = Begin(); it1 != End(); ++it1) {
                        for (auto it2 = Begin(); it2 != End(); ++it2) {
                            if (it1 != it2 && *it1 == *it2) {
                                return true;
                            }
                        }
                    }
                    return false;
                }
                
                //  Void functions
                void Fill(const T &value) {
                    std::fill(Begin(), End(), value);
                }
                
                void Zero() {
                    Fill(0);
                }
                
                void Ones() {
                    Fill(1);
                }
                
                void Identity() {
                    if (!IsSquared()) {
                        throw std::logic_error("Unable to create an identity matrix in a non squere matrix.");
                    }
                    
                    Zero();
                    auto it_end = End() + m;
                    const size_t step = m + 1;
                    for (auto it = Begin(); it != it_end; it += step) {
                        *it = 1;
                    }
                }
                
                void Write(const std::string &filename, const std::string &path = "/tmp") const {
                    const std::string file_path(path + "/" + filename);
                    std::ofstream output(file_path.data());
                    output << *this;
                }
                
                Matrix<T> operator+(const Matrix<T> &matrix) const {
                    if (this->n != matrix.n || this->m != matrix.m) {
                        throw std::logic_error("Matrices must be of the same dimensions.");
                    }
                    
                    Matrix<T> new_matrix(n, m);
                    auto it_new = new_matrix.Begin();
                    auto it_matrix = matrix.Begin();
                    
                    for (auto it_this = this->Begin(); it_this != this->End(); ++it_this, ++it_matrix, ++it_new) {
                        *it_new = *it_this + *it_matrix;
                    }
                    
                    return new_matrix;
                }
                
                Matrix<T>& operator+=(const Matrix<T> &matrix) {
                    if (this->n != matrix.n || this->m != matrix.m) {
                        throw std::logic_error("Matrices must be of the same dimensions.");
                    }
                    
                    auto it_matrix = matrix.Begin();
                    for (auto it_this = this->Begin(); it_this != this->End(); ++it_this, ++it_matrix) {
                        *it_this += *it_matrix;
                    }
                    
                    return *this;
                }
                
                Matrix<T> operator-(const Matrix<T> &matrix) const {
                    if (this->n != matrix.n || this->m != matrix.m) {
                        throw std::logic_error("Matrices must be of the same dimensions.");
                    }
                    
                    Matrix<T> new_matrix(n, m);
                    auto it_new = new_matrix.Begin();
                    auto it_matrix = matrix.Begin();
                    
                    for (auto it_this = this->Begin(); it_this != this->End(); ++it_this, ++it_matrix, ++it_new) {
                        *it_new = *it_this - *it_matrix;
                    }
                    
                    return new_matrix;
                }
                
                Matrix<T>& operator-=(const Matrix<T> &matrix) {
                    if (this->n != matrix.n || this->m != matrix.m) {
                        throw std::logic_error("Matrices must be of the same dimensions.");
                    }
                    
                    auto it_matrix = matrix.Begin();
                    for (auto it_this = this->Begin(); it_this != this->End(); ++it_this, ++it_matrix) {
                        *it_this -= *it_matrix;
                    }
                    
                    return *this;
                }
                
                template<typename T2>
                Matrix<T> operator*(const T2 &value) const {
                    Matrix<T> new_matrix(n, m);
                    auto it_new = new_matrix.Begin();
                    
                    for (auto it_this = this->Begin(); it_this != this->End(); ++it_this, ++it_new) {
                        *it_new = *it_this * static_cast<T>(value);
                    }
                    
                    return new_matrix;
                }
                
                template<typename T2>
                Matrix<T> operator*(const Matrix<T2> &matrix) const {
                    if (this->m != matrix.n) {
                        throw std::logic_error("Matrices dimensions are not compatible.");
                    }
                    
                    Matrix<T> new_matrix(this->n, matrix.m, 0);
                    
                    for (size_t this_row = 0; this_row < this->n; ++this_row) {
                        for (size_t other_column = 0; other_column < matrix.m; ++other_column) {
                            for (size_t element = 0; element < this->m; ++element) {
                                new_matrix[this_row][other_column] += this->operator[](this_row)[element] * static_cast<T>(matrix[element][other_column]);
                            }
                        }
                    }
                    
                    return new_matrix;
                }
                
                template<typename T2>
                Matrix<T> &operator*=(const Matrix<T2> &matrix) {
                    *this = *this * matrix;
                    return *this;
                }
                
                template<typename T2>
                Matrix<T> operator*(const Vector<T2> &vector) {
                    if (m != 1) {
                        throw std::logic_error("Matrix and Vector are not compatible.");
                    }
                    
                    Matrix<T> new_matrix(n, vector.Size());
                    auto it_new = new_matrix.Begin();
                    
                    for (auto it_this = this->Begin(); it_this != this->End(); ++it_this) {
                        for (auto it_vector = vector.Begin(); it_vector != vector.End(); ++it_vector) {
                            *it_new = *it_this * *it_vector;
                            ++it_new;
                        }
                    }
                    
                    return new_matrix;
                }
                
                Matrix<T>& operator*=(const T &value) {
                    for (auto it = Begin(); it != End(); ++it) {
                        *it *= value;
                    }
                    return *this;
                }
                
                Matrix<T> operator/(const T &value) const {
                    Matrix<T> new_matrix(this->n, this->m);
                    auto it_new = new_matrix.Begin();
                    
                    for (auto it_this = this->Begin(); it_this != this->End(); ++it_this, ++it_new) {
                        *it_new = *it_this / value;
                    }
                    
                    return new_matrix;
                }
                
                Matrix<T>& operator/=(const T &value) {
                    for (auto it = Begin(); it != End(); ++it) {
                        *it /= value;
                    }
                    return *this;
                }
                
                Matrix<T> operator-() const {
                    Matrix<T> new_matrix(this->n, this->m);
                    std::transform(this->Begin(), this->End(), new_matrix.Begin(),
                                   [](const T &value) {
                                       return -value;
                                   });
                    
                    return new_matrix;
                }
                
                Matrix<T> Pow(const ssize_t &power) const {
                    if (!this->IsSquared()) {
                        throw std::logic_error("Matrix must be an square matrix");
                    }
                    
                    Matrix<T> new_matrix;
                    
                    switch (power) {
                        case -1:
                            throw std::logic_error("To compute the inverse matrix consider using LU::InverseMatrix");
                            break;
                            
                        case 0:
                            new_matrix.Resize(this->n, this->m);
                            new_matrix.Identity();
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
                static Matrix<T> Zero(const size_t &rows, const size_t &columns) {
                    return Matrix<T>(rows, columns, 0);
                }
                
                static Matrix<T> Ones(const size_t &rows, const size_t &columns) {
                    return Matrix<T>(rows, columns, 1);
                }
                
                static Matrix<T> Identity(const size_t &rows, const size_t &columns) {
                    Matrix<T> identity(rows, columns);
                    identity.Identity();
                    return identity;
                }
                
            };
            
            //  MARK: - Extra funcions
            template <typename T>
            Matrix<T> Transpose(const Vector<T> &vector) {
                Matrix<T> matrix(vector.Size(), 1);
                std::copy(vector.Begin(), vector.End(), matrix.Begin());
                return matrix;
            }
        
        } /* namespace math */
    } /* namespace containers */
} /* namespace cda */

//  MARK: - Extra operators
template <typename T>
cda::math::containers::Matrix<T> operator*(const T &value,
                                           const cda::math::containers::Matrix<T> &matrix) {
    
    cda::math::containers::Matrix<T> tmp(matrix.Rows(), matrix.Columns());
    auto it_tmp = tmp.Begin();
    
    for (auto it_matrix = matrix.Begin(); it_matrix != matrix.End(); ++it_matrix, ++it_tmp) {
        *it_tmp = *it_matrix * value;
    }
    
    return tmp;
}

template <typename T>
cda::math::containers::Vector<T> operator*(const cda::math::containers::Vector<T> &vector,
                                           const cda::math::containers::Matrix<T> &matrix) {
    
    const auto &rows = matrix.Rows();
    if (rows != vector.Size()) {
        throw std::logic_error("The vector and the matrix are incompatible");
    }
    
    const auto &columns = matrix.Columns();
    cda::math::containers::Vector<T> new_vector(columns, 0);
    
    auto it_new = new_vector.Begin();
    for (size_t column = 0; column < columns; ++column, ++it_new) {
        auto it_vector = vector.Begin();
        for (size_t row = 0; row < rows; ++row, ++it_vector) {
            *it_new += *it_vector * matrix[row][column];
        }
    }
    
    return new_vector;
}

template <typename T>
cda::math::containers::Matrix<T> operator&&(const cda::math::containers::Matrix<T> &left_matrix,
                                            const cda::math::containers::Matrix<T> &right_matrix) {
    
    if (left_matrix.Rows() != right_matrix.Rows()) {
        throw std::logic_error("Both matrices must have the same number of rows");
    }
    
    cda::math::containers::Matrix<T> new_matrix(left_matrix.Rows(), left_matrix.Columns() + right_matrix.Columns());
    new_matrix.SetMatrix(0, 0, left_matrix);
    new_matrix.SetMatrix(0, left_matrix.Columns(), right_matrix);
    
    return new_matrix;
}

template <typename T>
void operator||(cda::math::containers::Matrix<T> &left_matrix,
                cda::math::containers::Matrix<T> &right_matrix) {
    
    const size_t left_matrix_columns = left_matrix.Columns() / 2;
    
    right_matrix = left_matrix.GetMatrix(0, left_matrix_columns, left_matrix.Rows(), left_matrix.Columns() - left_matrix_columns);
    left_matrix.Resize(left_matrix.Rows(), left_matrix_columns);
}

template <typename T>
void operator>>(std::istream &input,
                cda::math::containers::Matrix<T> &matrix) {
    
    if (!input) {
        throw std::logic_error("Input is not avaiable");
    }
    
    size_t initial_size = 100;
    matrix.Resize(initial_size, 1);
    
    size_t rows = 0, columns = 0;
    
    bool is_file_first_row = true;
    char comma;
    std::string line, cell;
    size_t element = 0;
    while (std::getline(input, line)) {
        
        std::stringstream line_stream(line);
        while (line_stream >> matrix[element][0]) {
            
            ++element;
            
            if (element >= matrix.Rows()) {
                initial_size *= 2;
                matrix.Resize(matrix.Rows() + initial_size, 1);
            }
            
            if (is_file_first_row) {
                ++columns;
            }
            
            line_stream >> comma;
        }
        
        if (is_file_first_row) {
            is_file_first_row = false;
        }
        
        ++rows;
    }
    
    if (rows * columns != element) {
        throw std::out_of_range("The size of the matrix could not be determined properly.");
    }
    
    matrix.Resize(element, 1);
    matrix.ChangeDimensions(rows, columns);
}

template <typename T>
std::ostream& operator<<(std::ostream &output,
                         const cda::math::containers::Matrix<T> &matrix) {
    
    if (output.rdbuf() == std::cout.rdbuf()) {
    //        out.width();
    //        out << fixed;
    //        out.fill(' ');
    //        out.precision(6);
    //        if (M.Rows() == 1 && M.Columns() == 1) {
    //            out << M(0) << endl;
    //        } else if (M.Rows() == 1) {
    //            out << "[" << M(0);
    //            for (int j=1; j<M.Columns(); j++) {
    //                out << "\t" << M(j);
    //            }
    //            out << "]\n";
    //        } else {
    //            out << "⎡" << M(0);
    //            for (int j=1; j<M.Columns(); j++) {
    //                out << right << "\t" << M(j);
    //            }
    //            out.width();
    //            out << "\t⎤\n";
    //            for (int i=1; i<M.Rows()-1; i++) {
    //                out << "⎢" << M(i*M.Columns());
    //                for (int j=1; j<M.Columns(); j++) {
    //                    out << "\t" << M(i*M.Columns() + j);
    //                }
    //                out << "\t⎥\n";
    //            }
    //            out << "⎣" << M((M.Rows()-1)*M.Columns());
    //            for (int j=1; j<M.Columns(); j++) {
    //                out << "\t" << M((M.Rows()-1)*M.Columns() + j);
    //            }
    //            out << "\t⎦\n";
    //        }
    //        out.precision();
        } else {
            for (size_t row = 0; row < matrix.Rows(); ++row) {
                for (size_t column = 0; column < matrix.Columns(); ++column) {
                    std::cout << matrix[row][column] << ";";
                }
                std::cout << std::endl;
            }
        }
    
    return output;
}
