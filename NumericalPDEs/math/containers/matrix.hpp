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
#include <stdexcept>
#include <type_traits>

#include "../algorithms/find.hpp"
#include "../algorithms/factorization/lu.hpp"


namespace cda {
    namespace math {
        namespace containers {
        
            //  --- MATRIX CLASS ---
            template <typename T>
            class Matrix {
            private:
                
                typedef typename std::conditional<std::is_floating_point<T>::value, T, float>::type ValueTypeLU;
                
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
                
                typedef T ValueType;
                
                Matrix(const size_t &rows = 0, const size_t &columns = 0) :
                n(rows), m(columns), size(rows * columns),
                a(nullptr), it_end(nullptr) {
                    AllocateMemory(size);
                }
                
                Matrix(const size_t &rows, const size_t &columns, const ValueType &value) :
                n(rows), m(columns), size(rows * columns),
                a(nullptr), it_end(nullptr) {
                    AllocateMemory(size);
                    std::fill(this->Begin(), this->End(), value);
                }
                
                Matrix(const Matrix<ValueType> &matrix) :
                n(matrix.n), m(matrix.m), size(matrix.size),
                a(nullptr), it_end(nullptr) {
                    AllocateMemory(size);
                    std::copy(matrix.Begin(), matrix.End(), this->Begin());
                }
                
                template<typename ValueType2>
                Matrix(const Matrix<ValueType2> &matrix) :
                n(matrix.Rows()), m(matrix.Columns()), size(matrix.Size()),
                a(nullptr), it_end(nullptr) {
                    AllocateMemory(size);
                    std::copy(matrix.Begin(), matrix.End(), this->Begin());
                }
                
                Matrix(Matrix<ValueType> &&matrix) :
                n(matrix.n), m(matrix.m), size(matrix.size),
                a(matrix.a), it_end(matrix.it_end) {
                    matrix.n = matrix.m = matrix.size = 0;
                    matrix.a = matrix.it_end = nullptr;
                }
                
                template <size_t size>
                Matrix(const size_t &rows, const size_t &columns, const ValueType (&values)[size]):
                n(rows), m(columns), size(rows * columns),
                a(nullptr), it_end(nullptr) {
                    if (this->size != size) {
                        throw std::logic_error("Sizes do not match");
                    }
                    AllocateMemory(size);
                    std::copy(values, values + size, this->Begin());
                }
                
                template <size_t rows, size_t columns>
                Matrix(const ValueType (&values)[rows][columns]):
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
                
                ValueType *Begin() const {
                    return a;
                }
                
                ValueType *End() const {
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
                            std::fill_n(a + size, new_size - size, static_cast<ValueType>(0));
                        }
                        
                        n = rows;
                        size = new_size;
                    } else {
                        
                        Matrix<ValueType> tmp(rows, columns);
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
                
                Matrix<ValueType> &operator=(const Matrix<ValueType> &matrix) {
                    if (this != &matrix) {
                        Resize(matrix.n, matrix.m);
                        std::copy(matrix.Begin(), matrix.End(), this->Begin());
                    }
                    
                    return *this;
                }
                
                Matrix<ValueType> &operator=(Matrix<ValueType> &&matrix) {
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
                
                bool operator==(const Matrix<ValueType> &matrix) const {
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
                
                bool operator!=(const Matrix<ValueType> &matrix) const {
                    return !this->operator==(matrix);
                }
                
                Matrix<ValueType> GetRow(const size_t &row) const {
                    if (row >= n) {
                        throw std::out_of_range("Index out of bounds");
                    }
                    
                    Matrix<ValueType> tmp(1, m);
                    
                    const auto &it_row = this->operator[](row);
                    std::copy(it_row, it_row + m, tmp.Begin());
                    
                    return tmp;
                }
                
                Matrix<ValueType> GetColumn(const size_t &column) const {
                    if (column >= m) {
                        throw std::out_of_range("Index out of bounds");
                    }
                    
                    Matrix<ValueType> tmp(n, 1);
                    auto it_tmp = tmp.Begin();
                    
                    for (size_t i = 0; i < n; ++i, ++it_tmp) {
                        *it_tmp = this->operator[](i)[column];
                    }
                    
                    return tmp;
                }
                
                Matrix<ValueType> GetMatrix(const size_t &row, const size_t &column,
                                    const size_t &number_of_rows, const size_t &number_of_columns) const {
                    if (row + number_of_rows > n || column + number_of_columns > m) {
                        throw std::out_of_range("Index out of bounds");
                    }
                    
                    Matrix<ValueType> tmp(number_of_rows, number_of_columns);
                    
                    auto it_this = Begin() + row * m + column;
                    for (auto it_tmp = tmp.Begin(); it_tmp != tmp.End(); it_tmp += number_of_columns, it_this += m) {
                        std::copy(it_this, it_this + number_of_columns, it_tmp);
                    }
                    
                    return tmp;
                }
                
                Matrix<ValueType> GetMatrix(const size_t &row, const size_t &column) const {
                    return this->GetMatrix(row, column, n - row, m - column);
                }
                
                Vector<ValueType> GetDiagonal() const {
                    if (!this->IsSquare()) {
                        throw std::logic_error("Matrix must be an square matrix");
                    }
                    
                    Vector<ValueType> diagonal(n);
                    
                    auto it_this = Begin();
                    for (auto it_diagonal = diagonal.Begin(); it_diagonal != diagonal.End(); ++it_diagonal, it_this += m + 1) {
                        *it_diagonal = *it_this;
                    }
                    
                    return diagonal;
                }
                
                Vector<ValueType> GetRowAsVector(const size_t &row, const size_t &from_column = 0) const {
                    if (row >= this->n || from_column >= this->m) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    Vector<ValueType> vector(this->m - from_column);
                    vector.Copy(vector.Size(), this->operator[](row) + from_column);
                    
                    return vector;
                }
                
                Vector<ValueType> GetColumnAsVector(const size_t &column, const size_t &from_row = 0) const {
                    if (column >= this->m || from_row >= this->n) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    Vector<ValueType> vector(this->n - from_row);
                    
                    auto it_this = this->Begin() + this->m * from_row + column;
                    for (auto it_vector = vector.Begin(); it_vector != vector.End(); ++it_vector, it_this += m) {
                        *it_vector = *it_this;
                    }
                    
                    return vector;
                }
                
                //  Sets
                void SetColumn(const size_t &column, const Matrix<ValueType> &matrix) {
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
                
                void SetColumn(const size_t &column, const Vector<ValueType> &vector) {
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
                
                void SetRow(const size_t &row, const Matrix<ValueType> &matrix) {
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
                
                void SetRow(const size_t &row, const Vector<ValueType> &vector) {
                    if (row >= this->n) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    if (vector.Size() != this->m) {
                        throw std::logic_error("The vector and the row of the matrix must have the same number of elements");
                    }
                    
                    auto it_this = this->operator[](row);
                    std::copy(vector.Begin(), vector.End(), it_this);
                }
                
                void SetMatrix(const size_t &row, const size_t &column, const Matrix<ValueType> &matrix) {
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
                
                Matrix<ValueType> SumRows() const {
                    Matrix<ValueType> sum_rows(n, 1, 0);
                    
                    auto it_this = this->Begin();
                    for (auto it_sum = sum_rows.Begin(); it_sum != sum_rows.End(); ++it_sum) {
                        auto it_this_row_end = it_this + this->m;
                        for (; it_this != it_this_row_end; ++it_this) {
                            *it_sum += *it_this;
                        }
                    }
                    
                    return sum_rows;
                }
                
                Vector<ValueType> SumRowsAsVector() const {
                    return SumRows().GetColumnAsVector(0);
                }
                
                ValueType SumRow(const size_t &row) const {
                    if (row >= this->n) {
                        throw std::out_of_range("Index out of range");
                    }
                    
                    ValueType sum = 0;
                    
                    const auto &it_end = this->operator[](row) + this->m;
                    for (auto it_row = this->operator[](row); it_row != it_end; ++it_row) {
                        sum += *it_row;
                    }
                    
                    return sum;
                }
                
                Matrix<ValueType> SumColumns() const {
                    Matrix<ValueType> sum_columns(1, m, 0);
                    
                    for (auto it_this = this->Begin(); it_this != this->End();) {
                        for (auto it_sum = sum_columns.Begin(); it_sum != sum_columns.End(); ++it_sum, ++it_this) {
                            *it_sum += *it_this;
                        }
                    }
                    
                    return sum_columns;
                }
                
                Vector<ValueType> SumColumnsAsVector() const {
                    return SumColumns().GetRowAsVector(0);
                }
                
                ValueType SumColumn(const size_t &column) const {
                    if (column >= this->m) {
                        throw std::out_of_range("Index out of range");
                    }
                    
                    ValueType sum = 0;
                    
                    const auto &it_column_end = this->End() + column;
                    for (auto it_column = this->Begin() + column; it_column != it_column_end; it_column += this->m) {
                        sum += *it_column;
                    }
                    
                    return sum;
                }
                
                ValueType MaximumElement() const {
                    return algorithms::find::MaximumElement(Begin(), End());
                }
                
                ValueType AbsoluteMaximumElement() const {
                    return algorithms::find::AbsoluteMaximumElement(Begin(), End());
                }
                
                ValueType AbsoluteMaximumElementWithSign() const {
                    return algorithms::find::AbsoluteMaximumElementWithSign(Begin(), End());
                }
                
                ValueType MinimumElement() const {
                    return algorithms::find::MinimumElement(Begin(), End());
                }
                
                ValueType AbsoluteMinimumElement() const {
                    return algorithms::find::AbsoluteMinimumElement(Begin(), End());
                }
                
                ValueType AbsoluteMinimumElementWithSign() const {
                    return algorithms::find::AbsoluteMinimumElementWithSign(Begin(), End());
                }
                
                const ValueType &At(const size_t &row, const size_t &column) const {
                    if (row >= n || column >= m) {
                        throw std::out_of_range("Indexes out of range");
                    }
                    return this->operator[](row)[column];
                }
                
                const ValueType *operator[](const size_t &row) const {
                    return &a[row * m];
                }
                
                ValueType &At(const size_t &row, const size_t &column) {
                    if (row >= n || column >= m) {
                        throw std::out_of_range("Indexes out of range");
                    }
                    return this->operator[](row)[column];
                }
                
                ValueType *operator[](const size_t &row) {
                    return &a[row * m];
                }
                
                Matrix<ValueType> Transpose() const {
                    Matrix<ValueType> new_matrix(this->m, this->n);
                    
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
                
                bool IsSquare() const {
                    return n == m;
                }
                
                bool HasDuplicate(const ValueType &accuracy) const {
                    ValueType distance;
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
                void Fill(const ValueType &value) {
                    std::fill(Begin(), End(), value);
                }
                
                void Zero() {
                    Fill(0);
                }
                
                void Ones() {
                    Fill(1);
                }
                
                void Identity() {
                    if (!IsSquare()) {
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
                
                Matrix<ValueType> operator+(const Matrix<ValueType> &matrix) const {
                    if (this->n != matrix.n || this->m != matrix.m) {
                        throw std::logic_error("Matrices must be of the same dimensions.");
                    }
                    
                    Matrix<ValueType> new_matrix(n, m);
                    auto it_new = new_matrix.Begin();
                    auto it_matrix = matrix.Begin();
                    
                    for (auto it_this = this->Begin(); it_this != this->End(); ++it_this, ++it_matrix, ++it_new) {
                        *it_new = *it_this + *it_matrix;
                    }
                    
                    return new_matrix;
                }
                
                Matrix<ValueType>& operator+=(const Matrix<ValueType> &matrix) {
                    if (this->n != matrix.n || this->m != matrix.m) {
                        throw std::logic_error("Matrices must be of the same dimensions.");
                    }
                    
                    auto it_matrix = matrix.Begin();
                    for (auto it_this = this->Begin(); it_this != this->End(); ++it_this, ++it_matrix) {
                        *it_this += *it_matrix;
                    }
                    
                    return *this;
                }
                
                Matrix<ValueType> operator-(const Matrix<ValueType> &matrix) const {
                    if (this->n != matrix.n || this->m != matrix.m) {
                        throw std::logic_error("Matrices must be of the same dimensions.");
                    }
                    
                    Matrix<ValueType> new_matrix(n, m);
                    auto it_new = new_matrix.Begin();
                    auto it_matrix = matrix.Begin();
                    
                    for (auto it_this = this->Begin(); it_this != this->End(); ++it_this, ++it_matrix, ++it_new) {
                        *it_new = *it_this - *it_matrix;
                    }
                    
                    return new_matrix;
                }
                
                Matrix<ValueType>& operator-=(const Matrix<ValueType> &matrix) {
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
                Matrix<ValueType> operator*(const T2 &value) const {
                    Matrix<ValueType> new_matrix(n, m);
                    auto it_new = new_matrix.Begin();
                    
                    for (auto it_this = this->Begin(); it_this != this->End(); ++it_this, ++it_new) {
                        *it_new = *it_this * static_cast<ValueType>(value);
                    }
                    
                    return new_matrix;
                }
                
                template<typename T2>
                Matrix<ValueType> operator*(const Matrix<T2> &matrix) const {
                    if (this->m != matrix.n) {
                        throw std::logic_error("Matrices dimensions are not compatible.");
                    }
                    
                    Matrix<ValueType> new_matrix(this->n, matrix.m, 0);
                    
                    for (size_t this_row = 0; this_row < this->n; ++this_row) {
                        for (size_t other_column = 0; other_column < matrix.m; ++other_column) {
                            for (size_t element = 0; element < this->m; ++element) {
                                new_matrix[this_row][other_column] += this->operator[](this_row)[element] * static_cast<ValueType>(matrix[element][other_column]);
                            }
                        }
                    }
                    
                    return new_matrix;
                }
                
                template<typename T2>
                Matrix<ValueType> &operator*=(const Matrix<T2> &matrix) {
                    *this = *this * matrix;
                    return *this;
                }
                
                template<typename T2>
                Matrix<ValueType> operator*(const Vector<T2> &vector) {
                    if (m != 1) {
                        throw std::logic_error("Matrix and Vector are not compatible.");
                    }
                    
                    Matrix<ValueType> new_matrix(n, vector.Size());
                    auto it_new = new_matrix.Begin();
                    
                    for (auto it_this = this->Begin(); it_this != this->End(); ++it_this) {
                        for (auto it_vector = vector.Begin(); it_vector != vector.End(); ++it_vector) {
                            *it_new = *it_this * *it_vector;
                            ++it_new;
                        }
                    }
                    
                    return new_matrix;
                }
                
                Matrix<ValueType>& operator*=(const ValueType &value) {
                    for (auto it = Begin(); it != End(); ++it) {
                        *it *= value;
                    }
                    return *this;
                }
                
                Matrix<ValueType> operator/(const ValueType &value) const {
                    Matrix<ValueType> new_matrix(this->n, this->m);
                    auto it_new = new_matrix.Begin();
                    
                    for (auto it_this = this->Begin(); it_this != this->End(); ++it_this, ++it_new) {
                        *it_new = *it_this / value;
                    }
                    
                    return new_matrix;
                }
                
                Matrix<ValueType>& operator/=(const ValueType &value) {
                    for (auto it = Begin(); it != End(); ++it) {
                        *it /= value;
                    }
                    return *this;
                }
                
                Matrix<ValueType> operator-() const {
                    Matrix<ValueType> new_matrix(this->n, this->m);
                    std::transform(this->Begin(), this->End(), new_matrix.Begin(),
                                   [](const ValueType &value) {
                                       return -value;
                                   });
                    
                    return new_matrix;
                }
                
                ValueType Determinant() const {
                    return algorithms::factorization::LU<Matrix, ValueTypeLU>::Determinant(*this);
                }
                
                Matrix<ValueType> Pow(const ssize_t &power) const {
                    if (!this->IsSquare()) {
                        throw std::logic_error("Matrix must be an square matrix");
                    }
                    
                    Matrix<ValueType> new_matrix;
                    
                    switch (power) {
                        case -1:
                            new_matrix = algorithms::factorization::LU<Matrix, ValueTypeLU>::InverseMatrix(*this);
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
                static Matrix<ValueType> Zero(const size_t &rows, const size_t &columns) {
                    return Matrix<ValueType>(rows, columns, 0);
                }
                
                static Matrix<ValueType> Ones(const size_t &rows, const size_t &columns) {
                    return Matrix<ValueType>(rows, columns, 1);
                }
                
                static Matrix<ValueType> Identity(const size_t &rows) {
                    Matrix<ValueType> identity(rows, rows);
                    identity.Identity();
                    return identity;
                }
                
            };
            
            //  MARK: - Extra funcions
            template <typename ValueType>
            Matrix<ValueType> Transpose(const Vector<ValueType> &vector) {
                Matrix<ValueType> matrix(vector.Size(), 1);
                std::copy(vector.Begin(), vector.End(), matrix.Begin());
                return matrix;
            }
        
        } /* namespace math */
    } /* namespace containers */
} /* namespace cda */

//  MARK: - Extra operators
template <typename ValueType>
cda::math::containers::Matrix<ValueType> operator*(const ValueType &value,
                                           const cda::math::containers::Matrix<ValueType> &matrix) {
    
    cda::math::containers::Matrix<ValueType> tmp(matrix.Rows(), matrix.Columns());
    auto it_tmp = tmp.Begin();
    
    for (auto it_matrix = matrix.Begin(); it_matrix != matrix.End(); ++it_matrix, ++it_tmp) {
        *it_tmp = *it_matrix * value;
    }
    
    return tmp;
}

template <typename ValueType>
cda::math::containers::Vector<ValueType> operator*(const cda::math::containers::Vector<ValueType> &vector,
                                           const cda::math::containers::Matrix<ValueType> &matrix) {
    
    const auto &rows = matrix.Rows();
    if (rows != vector.Size()) {
        throw std::logic_error("The vector and the matrix are incompatible");
    }
    
    const auto &columns = matrix.Columns();
    cda::math::containers::Vector<ValueType> new_vector(columns, 0);
    
    auto it_new = new_vector.Begin();
    for (size_t column = 0; column < columns; ++column, ++it_new) {
        auto it_vector = vector.Begin();
        for (size_t row = 0; row < rows; ++row, ++it_vector) {
            *it_new += *it_vector * matrix[row][column];
        }
    }
    
    return new_vector;
}

template <typename ValueType>
cda::math::containers::Matrix<ValueType> operator&&(const cda::math::containers::Matrix<ValueType> &left_matrix,
                                            const cda::math::containers::Matrix<ValueType> &right_matrix) {
    
    if (left_matrix.Rows() != right_matrix.Rows()) {
        throw std::logic_error("Both matrices must have the same number of rows");
    }
    
    cda::math::containers::Matrix<ValueType> new_matrix(left_matrix.Rows(), left_matrix.Columns() + right_matrix.Columns());
    new_matrix.SetMatrix(0, 0, left_matrix);
    new_matrix.SetMatrix(0, left_matrix.Columns(), right_matrix);
    
    return new_matrix;
}

template <typename ValueType>
void operator||(cda::math::containers::Matrix<ValueType> &left_matrix,
                cda::math::containers::Matrix<ValueType> &right_matrix) {
    
    const size_t left_matrix_columns = left_matrix.Columns() / 2;
    
    right_matrix = left_matrix.GetMatrix(0, left_matrix_columns, left_matrix.Rows(), left_matrix.Columns() - left_matrix_columns);
    left_matrix.Resize(left_matrix.Rows(), left_matrix_columns);
}

template <typename ValueType>
void operator>>(std::istream &input,
                cda::math::containers::Matrix<ValueType> &matrix) {
    
    if (!input) {
        throw std::logic_error("Input is not avaiable");
    }
    
    size_t initial_size = 100;
    matrix.Resize(initial_size, 1);
    
    size_t rows = 0, columns = 0;
    
    bool is_first_row = true;
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
    
    matrix.Resize(element, 1);
    matrix.ChangeDimensions(rows, columns);
}

template <typename ValueType>
std::ostream& operator<<(std::ostream &output,
                         const cda::math::containers::Matrix<ValueType> &matrix) {
    
    const size_t rows = matrix.Rows();
    const size_t columns = matrix.Columns();
    
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
        for (size_t row = 0; row < rows; ++row) {
            for (size_t column = 0; column < columns; ++column) {
                output << matrix[row][column] << ";";
            }
            output << std::endl;
        }
    }

    return output;
}
