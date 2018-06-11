//
//  matrix.hpp
//  Vibration Normal Modes 2D
//
//  Created by Carlos David on 04/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#ifndef math_containers_matrix_hpp
#define math_containers_matrix_hpp

#include "vector.hpp"

#include <utility>

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
                    if (void *mem = std::realloc(a, size * sizeof(T))) {
                        a = static_cast<T *>(mem);
                        it_end = a + size;
                    } else {
                        throw std::bad_alloc();
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
                    
                    const auto &it_row = a + row * m;
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
                        *it_tmp = a[i*m + column];
                    }
                    
                    return tmp;
                }
                
                Matrix<T> GetMatrix(const size_t &row, const size_t &column, const size_t &number_of_rows, const size_t &number_of_columns) const {
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
                    if (from_column >= this->m) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    Vector<T> vector(this->m - from_column);
                    vector.Copy(vector.Size(), this->operator[](row) + from_column);
                    
                    return vector;
                }
                
                Vector<T> GetColumnAsVector(const size_t &column, const size_t &from_row = 0) const {
                    if (from_row >= this->n) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    Vector<T> vector(this->n - from_row);
                    
                    auto it_this = this->Begin() + this->m * from_row;
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
                    
                    auto it_this = this->Begin() + row * this->m;
                    std::copy(matrix.Begin(), matrix.End(), it_this);
                }
                
                void SetRow(const size_t &row, const Vector<T> &vector) {
                    if (row >= this->n) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    if (vector.Size() != this->m) {
                        throw std::logic_error("The vector and the row of the matrix must have the same number of elements");
                    }
                    
                    auto it_this = this->Begin() + row * this->m;
                    std::copy(vector.Begin(), vector.End(), it_this);
                }
                
                void SetMatrix(const size_t &row, const size_t &column, const Matrix<T> &matrix) {
                    if (row >= this->n || column >= this->m) {
                        throw std::out_of_range("Index out of bounds.");
                    }
                    
                    const auto number_of_rows = std::min(this->n - row, matrix.Rows());
                    const auto number_of_columns = std::min(this->m - column, matrix.Columns());
                    
                    auto it_matrix = matrix.Begin();
                    const auto it_this_end = this->operator[](row) + number_of_rows;
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
                
                T SumRow(const size_t &row) const {
                    T sum = 0;
                    
                    const auto &it_end = this->operator[](row) + this->m;
                    for (auto it_row = this->operator[](row); it_row != it_end; ++it_row) {
                        sum += *it_row;
                    }
                    
                    return sum;
                }
                
                T SumColumn(const size_t &column) const {
                    T sum = 0;
                    
                    const auto &it_column_end = this->End() + column;
                    for (auto it_column = this->Begin() + column; it_column != it_column_end; it_column += this->m) {
                        sum += *it_column;
                    }
                    
                    return sum;
                }
                
                T max() const;                                                          //  Devuelve el máximo de la matriz.
                T maxAbs() const;                                                       //  Devuelve el máximo absoluto de la matriz.
                T maxAbs_sig() const;                                                   //  Devuelve el máximo absoluto de la matriz, pero con su signo.
                T min() const;                                                          //  Devuelve el mínimo de la matriz.
                T det();                                                                //  Calcula el determinante de una matriz.
                
                const T &At(const size_t &row, const size_t &column) const {
                    if(row >= n || column >= m) {
                        throw std::out_of_range("Indexes out of range");
                    }
                    return a[row * m + column];
                }
                
                const T *operator[](const size_t &row) const {
                    return &a[row * m];
                }
                
                T &At(const size_t &row, const size_t &column) {
                    if (row >= n || column >= m) {
                        throw std::out_of_range("Indexes out of range");
                    }
                    return a[row * m + column];
                }
                
                T *operator[](const size_t &row) {
                    return &a[row * m];
                }
                
                //  Returns a Matrix / a vector
                Matrix<T> sumRows();                                                      //  Suma todos los elementos de las filas y las devuelve en una matriz columna.
                Matrix<T> sumColumns();                                                   //  Suma todos los elementos de las columnas y las devuelve en una matriz fila.
                Matrix<T> transpose();                                                    //  Transpone la matriz.
                
                Vector<T> sumRowsV();                                                     //  Suma todos los elementos de las filas y las devuelve en un vector.
                Vector<T> sumColumnsV();                                                  //  Suma todos los elementos de las columnas y las devuelve en un vector.
                
                //  Returns a bool
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
                
                bool duplicate(T range);                                                //  Comprueba si hay elementos duplicados (próximos);
                bool duplicate();                                                       //  Comprueba si hay elementos duplicados.
                
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
                
                void write(const std::string& path, const std::string& filename);                 //  Escribe la matriz en un fichero.
                void write(const std::string& filename);
                
                
                
                //  --- MÉTODO LU ---
                Matrix<T> LU();                                                           //  Descompone una matriz en dos por el método LU.
                Matrix<T> U();                                                            //  Devuelve la matriz U del método LU.
                Matrix<T> L();                                                            //  Devuelve la matriz L del método LU.
                Vector<T> solveLU(const Vector<T>& B);                                      //  Resuelve un sistema de ecuaciones por el método LU.
                Vector<T> solveLU3d(const Vector<T>& B);                                    //  Resuelve un sistema de ecuaciones tridiagonal por el método LU.
                Vector<T> solveGS3d(const Vector<T>& B, T err);                             //  Resuelve un sistema de ecuaciones tridiagonal por el método Gauss-Seidel.
                
                //  --- MÉTODO QR ---
                Matrix<T> QR(unsigned char QorR);                                         //  Descompone una matriz en dos por el método QR.
                Matrix<T> eigenVectors(int maxIte, unsigned char opt);                    //  Calcula los autovectores a partir de los autovalores (QR).
                Matrix<T> eigenVectors(unsigned char opt);                                //  Calcula los autovectores a partir de los autovaleres (QR).
                Matrix<T> eigenVectors();                                                 //  Calcula los autovectores a partir de los autovalores (QR).
                Vector<T> eigenVector(T eigenValue, int maxIte, unsigned char opt);       //  Calcula el autovector correspondiente al autovector introducido.
                Vector<T> eigenVector(T eigenValue, unsigned char opt);                   //  Calcula el autovector correspondiente al autovector introducido.
                Vector<T> eigenVector(T eigenValue);                                      //  Calcula el autovector correspondiente al autovector introducido.
                Vector<T> eigenValues(int maxIte, T factor, unsigned char opt);           //  Calcula los autovalores de una matriz por el método QR.
                Vector<T> eigenValues(int maxIte, unsigned char opt);                     //  Calcula los autovalores de una matriz por el método QR.
                Vector<T> eigenValues(unsigned char opt);                                 //  Calcula los autovalores de una matriz por el método QR.
                Vector<T> eigenValues();                                                  //  Calcula los autovalores de una matraz por el método QR.
                
                
                
                //  --- OPERATORS ---
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
                Matrix<T> operator*(const Matrix<T2> &matrix) {
                    if (this->n != matrix.m || this->m != matrix.n) {
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
                Matrix<T> operator*(const Vector<T2> &vector) {
                    if (m != 1) {
                        throw std::logic_error("Matrix and Vector are not compatible.");
                    }
                    
                    Matrix<T> new_matrix(n, vector.Size());
                    auto it_new = new_matrix.Begin();
                    auto it_vector = vector.Begin();
                    
                    for (auto it_this = this->Begin(); it_this != this->End(); ++it_this, ++it_vector, ++it_new) {
                        *it_new = *it_this * *it_vector;
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
                
                Matrix<T> operator-();                                                    //  Cambia el signo de todos los elementos de la matriz.
                Matrix<T>& operator*=(const Matrix<T>& M);                                  //  Producto de una matriz por otra sobre sí misma.
                Matrix<T> operator^(const int exp);                                       //  Eleva a una potencia la matriz ó calcula su inversa.
                
                //  --- STATIC METHODS ---
                static Matrix<T> Zero(const size_t &rows, const size_t &columns) {
                    return Matrix<T>(rows, columns, 0);
                }
                
                static Matrix<T> Ones(const size_t &rows, const size_t &columns) {
                    return Matrix<T>(rows, columns, 1);
                }
                
                static Matrix<T> Identity(const size_t &rows, const size_t &columns) {
                    if (rows != columns) {
                        throw std::logic_error("Unable to create an identity matrix in a non squere matrix.");
                    }
                    
                    Matrix<T> identity(rows, columns);
                    identity.Identity();
                    return identity;
                }
                
            };
            
            template <typename T>
            Matrix<T> setdiff(Matrix<T>& A, const Matrix<T>& B, const int& reps);                     //  Devuelve los valores de A que no se encuentran en B.
            template <typename T>
            Matrix<T> setdiff(Matrix<T>& A, const Matrix<T>& B);                                      //  Devuelve los valores de A que no se encuentran en B.
        
        } /* namespace math */
    } /* namespace containers */
} /* namespace cda */

//  --- MORE OPERATORS ---
template <typename T>
inline cda::math::containers::Matrix<T> operator*(const T &value, const cda::math::containers::Matrix<T> &matrix) {
    cda::math::containers::Matrix<T> tmp(matrix.Rows(), matrix.Columns());
    auto it_tmp = tmp.Begin();
    
    for (auto it_matrix = matrix.Begin(); it_matrix != matrix.End(); ++it_matrix, ++it_tmp) {
        *it_tmp = *it_matrix * value;
    }
    
    return tmp;
}

template<typename T>
inline cda::math::containers::Vector<T> operator*(const cda::math::containers::Vector<T> &vector,
                                                  const cda::math::containers::Matrix<T> &matrix) {
    const auto &rows = matrix.Rows();
    if (rows != vector.Size()) {
        throw std::logic_error("The vector and the matrix are incompatible");
    }
    
    const auto &columns = matrix.Columns();
    cda::math::containers::Vector<T> new_vector(columns, 0);
    
    // TODO: Check this operation
    auto it_new = new_vector.Begin();
    for (size_t column = 0; column < columns; ++column, ++it_new) {
        auto it_vector = vector.Begin();
        for (size_t row = 0; row < rows; ++row, ++it_vector) {
            *it_new += *it_vector * matrix[row][column];
        }
    }
    
    return new_vector;
}

template <typename T> inline cda::math::containers::Matrix<T>
operator&&(const cda::math::containers::Matrix<T>& Mi, const cda::math::containers::Matrix<T>& Md);                           //  Devuelve una matriz que tendrá unidas dos matrices.
template <typename T> inline void
operator||(cda::math::containers::Matrix<T>& Mi, cda::math::containers::Matrix<T>& Md);                                       //  Separa una matriz en dos.
template <typename T> inline void
operator>>(std::istream& in, cda::math::containers::Matrix<T>& M);                                        //  Importa una matriz desde un fichero.
template <typename T> inline std::ostream&
operator<<(std::ostream& out, const cda::math::containers::Matrix<T>& M);                                 //  Exporta la matriz a un ostream.

#endif /* math_containers_matrix_hpp */
