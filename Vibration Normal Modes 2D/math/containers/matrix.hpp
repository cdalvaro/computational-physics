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

namespace cda {
    namespace math {
        namespace containers {
        
            //  --- MATRIX CLASS ---
            template <typename T>
            class Matrix {
            private:
                size_t n, m, size;
                T *a;
                
                void AllocateMemory(const size_t &size) {
                    if (void *mem = std::realloc(a, size * sizeof(T))) {
                        a = static_cast<T *>(mem);
                    } else {
                        throw std::bad_alloc();
                    }
                }
                
            public:
                
                Matrix(const size_t &rows = 0, const size_t &columns = 0) :
                n(rows), m(columns), size(rows * columns), a(nullptr) {
                    AllocateMemory(size);
                }
                
                Matrix(const size_t &rows, const size_t &columns, const T &value) :
                n(rows), m(columns), size(rows * columns), a(nullptr) {
                    AllocateMemory(size);
                    std::fill(a, a + size, value);
                }
                
                Matrix(const Matrix<T> &matrix) :
                n(matrix.n), m(matrix.m), size(matrix.size), a(nullptr) {
                    AllocateMemory(size);
                    std::copy(matrix.a, matrix.a + size, a);
                }
                
                Matrix(Matrix<T> &&matrix) :
                n(matrix.n), m(matrix.m), size(matrix.size), a(matrix.a) {
                    matrix.n = 0;
                    matrix.m = 0;
                    matrix.size = 0;
                    matrix.a = nullptr;
                }
                
                ~Matrix() {
                    if (a) {
                        std::free(a);
                        a = nullptr;
                        n = m = size = 0;
                    }
                }
                
                T *Begin() const {
                    return a;
                }
                
                T *End() const {
                    return a + size;
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
                            tmp.zero();
                        }
                        
                        if (size != 0) {
                            const auto &rows_to_copy = rows < n ? rows : n;
                            const auto &columns_to_copy = columns < m ? columns : m;
                            for (size_t row = 0; row < rows_to_copy; ++row) {
                                for (size_t column = 0; column < columns_to_copy; ++column) {
                                    tmp(row, column) = this->operator()(row, column);
                                }
                            }
                        }
                        
                        this->operator=(std::move(tmp));
                    }
                    
                }
                
                Matrix<T> &operator=(const Matrix<T> &matrix) {
                    if (this != &matrix) {
                        Resize(matrix.n, matrix.m);
                        std::copy(matrix.a, matrix.a + size, a);
                    }
                    
                    return *this;
                }
                
                Matrix<T> &operator=(Matrix<T> &&matrix) {
                    if (this != &matrix) {
                        std::free(a);
                        a = matrix.a;
                        n = matrix.n;
                        m = matrix.m;
                        size = matrix.size;
                        
                        matrix.a = nullptr;
                        matrix.n = 0;
                        matrix.m = 0;
                        matrix.size = 0;
                    }
                    
                    return *this;
                }
                
                Matrix<T> GetRow(const size_t &row) const {
                    if (row >= n) {
                        throw std::out_of_range("Index out of bounds");
                    }
                    
                    Matrix<T> tmp(1, m);
                    
                    const auto &it_row = a + row * m;
                    std::copy(it_row, it_row + m, tmp.a);
                    
                    return tmp;
                }
                
                Matrix<T> GetColumn(const size_t &column) const {
                    if (column >= m) {
                        throw std::out_of_range("Index out of bounds");
                    }
                    
                    Matrix<T> tmp(n, 1);
                    auto it_tmp = tmp.a;
                    
                    for (size_t i = 0; i < n; ++i, ++it_tmp) {
                        *it_tmp = a[i*m + column];
                    }
                    
                    return tmp;
                }
                
                Matrix<T> getMatrix(int _i, int _j, int rows, int columns);               //  Devuelve la submatriz (i,j) con tamaño (n,m).
                Matrix<T> getMatrix(int _i, int _j);                                      //  Devuelve la submatriz (i,j).
                
                Vector<T> getDiagonal();                                                  //  Devuelve la diagonal de una matriz.
                Vector<T> getRowV(int _i, int _j);                                        //  Devuelve como vector la fila i desde la columna j.
                Vector<T> getRowV(int _i);                                                //  Devuelve como vector la fila i.
                Vector<T> getColumnV(int _i, int _j);                                     //  Devuelve como vector la columna j desde la fila i.
                Vector<T> getColumnV(int _j);                                             //  Devuelve como vector la columna j.
                
                //  Sets
                void setColumn(int _i, int _j, const Matrix<T>& Mc);                      //  Introduce una matriz columna en la columna j desde la fila i.
                void setColumn(int _j, const Matrix<T>& Mc);                              //  Introduce una matriz columna en la columna j.
                void setColumnV(int _i, int _j, const Vector<T>& V);                      //  Introduce un vector en la columna j desde la fila i.
                void setColumnV(int _j, const Vector<T>& V);                              //  Introduce un vector en la columna j.
                void setRow(int _i, int _j, const Matrix<T>& Mr);                         //  Introduce una matriz fila en la fila i desde la columna j.
                void setRow(int _i, const Matrix<T>& Mr);                                 //  Introduce una matriz fila en la fila i.
                void setRowV(int _i, int _j, const Vector<T>& V);                         //  Introduce un vector en la fila i desde la columna j.
                void setRowV(int _i, const Vector<T>& V);                                 //  Introduce un vector en la fila i.
                void setMatrix(int _i, int _j, const Matrix<T>& M);                       //  Introduce una matriz en la posición (i,j).
                
                
                
                //  --- FUNCTIONS ---
                //  Returns an integer / a double / a float
                int dims() const;                                                       //  Devuelve las dimensiones de la matriz.
                int rows() const;                                                       //  Devuelve el número de filas de la matriz.
                int columns() const;                                                    //  Devuelve el número de columnas de la matriz.
                int elements() const;                                                   //  Devuelve el número de elementos de la matriz.
                
                T sumRow(int _i) const;                                                 //  Suma todos los elementos de la fila i.
                T sumColumn(int _j) const;                                              //  Suma todos los elementos de la columna j.
                T max() const;                                                          //  Devuelve el máximo de la matriz.
                T maxAbs() const;                                                       //  Devuelve el máximo absoluto de la matriz.
                T maxAbs_sig() const;                                                   //  Devuelve el máximo absoluto de la matriz, pero con su signo.
                T min() const;                                                          //  Devuelve el mínimo de la matriz.
                T det();                                                                //  Calcula el determinante de una matriz.
                T& operator()(int _i, int _j) const;                                    //  Cambia el valor del elemento (i,j). (Notación matricial).
                
                const T *operator[](const size_t &row) const {
                    return &a[row * m];
                }
                
                T *operator[](const size_t &row) {
                    return &a[row * m];
                }
                
                //  Returns a Matrix / a vector
                Matrix<T> sumRows();                                                      //  Suma todos los elementos de las filas y las devuelve en una matriz columna.
                Matrix<T> sumColumns();                                                   //  Suma todos los elementos de las columnas y las devuelve en una matriz fila.
                Matrix<T> transpose();                                                    //  Transpone la matriz.
                Matrix<T> zero(int _n, int _m);                                           //  Devuelve una matriz de tamaño (n,m) llena de ceros.
                Matrix<T> ones(int _n, int _m);                                           //  Devuelve una matriz de tamaño (n,m) llena de unos.
                Matrix<T> identity(int _n, int _m);                                       //  Devuelve una matriz identidad de tamaño (n,m).
                
                Vector<T> sumRowsV();                                                     //  Suma todos los elementos de las filas y las devuelve en un vector.
                Vector<T> sumColumnsV();                                                  //  Suma todos los elementos de las columnas y las devuelve en un vector.
                
                //  Returns a bool
                bool null();                                                            //  Comprueba si la matriz es una matriz de ceros.
                bool duplicate(T range);                                                //  Comprueba si hay elementos duplicados (próximos);
                bool duplicate();                                                       //  Comprueba si hay elementos duplicados.
                
                //  Void functions
                void Fill(const T &value) {
                    std::fill(a, a + size, value);
                }
                
                void zero() {
                    Fill(0);
                }
                
                void ones() {
                    Fill(1);
                }
                
                void identity();                                                        //  Combierte a la matriz en la identidad.
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
                Matrix<T>& operator=(const T* Array);                                     //  Introduce los elementos de una matriz a través de un array.
                Matrix<T> operator+(const Matrix<T>& M);                                    //  Definición de suma de dos matrices.
                Matrix<T>& operator+=(const Matrix<T>& M);                                  //  Suma la matriz con otra sobre sí misma.
                Matrix<T> operator-(const Matrix<T>& M);                                    //  Definición de diferencia de dos matrices.
                Matrix<T>& operator-=(const Matrix<T>& M);                                  //  Diferencia de la matriz con otra sobre sí misma.
                Matrix<T> operator*(const T& D);                                          //  Producto de la matriz por una constante.
                Matrix<T>& operator*=(const T& D);                                        //  Producto de la matriz sobre sí misma por una constante.
                Matrix<T> operator/(const T& D);                                          //  Cociente de la matriz por una constante.
                Matrix<T>& operator/=(const T& D);                                        //  Cociente de la matriz sobre sí misma por una constante.
                Matrix<T> operator-();                                                    //  Cambia el signo de todos los elementos de la matriz.
                Matrix<T> operator*(const Matrix<T>& M);                                    //  Definición del producto de dos matrices.
                Matrix<T>& operator*=(const Matrix<T>& M);                                  //  Producto de una matriz por otra sobre sí misma.
                Matrix<T> operator*(const Vector<T>& V);                                    //  Producto de una matriz por un vector.
                Matrix<T> operator^(const int exp);                                       //  Eleva a una potencia la matriz ó calcula su inversa.
            };
            
            //  --- OTHER FUNCTIONS ---
            template <class T>
            Matrix<T> zero(int _n, int _m);                                                       //  Devuelve una matriz de tamaño (n,m) llena de ceros.
            template <class T>
            Matrix<T> ones(int _n, int _m);                                                       //  Devuelve una matriz de tamaño (n,m) llena de unos.
            template <class T>
            Matrix<T> identity(int _n, int _m);                                                   //  Devuelve una matriz identidad de tamaño (n,m).
            template <class T>
            Matrix<T> setdiff(Matrix<T>& A, const Matrix<T>& B, const int& reps);                     //  Devuelve los valores de A que no se encuentran en B.
            template <class T>
            Matrix<T> setdiff(Matrix<T>& A, const Matrix<T>& B);                                      //  Devuelve los valores de A que no se encuentran en B.
        
        } /* namespace math */
    } /* namespace containers */
} /* namespace cda */

//  --- MORE OPERATORS ---
template <class T>
inline cda::math::containers::Matrix<T> operator*(const T &value, const cda::math::containers::Matrix<T> &matrix) {
    cda::math::containers::Matrix<T> tmp(matrix.rows(), matrix.columns());
    auto it_tmp = tmp.Begin();
    
    for (auto it_matrix = matrix.Begin(); it_matrix != matrix.End(); ++it_matrix, ++it_tmp) {
        *it_tmp = *it_matrix * value;
    }
    
    return tmp;
}

template <class T> inline cda::math::containers::Matrix<T>
operator&&(const cda::math::containers::Matrix<T>& Mi, const cda::math::containers::Matrix<T>& Md);                           //  Devuelve una matriz que tendrá unidas dos matrices.
template <class T> inline void
operator||(cda::math::containers::Matrix<T>& Mi, cda::math::containers::Matrix<T>& Md);                                       //  Separa una matriz en dos.
template <class T> inline void
operator>>(std::istream& in, cda::math::containers::Matrix<T>& M);                                        //  Importa una matriz desde un fichero.
template <class T> inline std::ostream&
operator<<(std::ostream& out, const cda::math::containers::Matrix<T>& M);                                 //  Exporta la matriz a un ostream.

#endif /* math_containers_matrix_hpp */
