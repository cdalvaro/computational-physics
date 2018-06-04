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
                int n, m, size;
                T *a;
                
            public:
                //  --- DEFINITION ---
                Matrix(int _n, int _m) :
                n(_n), m(_m), size(_n * _m), a(new T[size]) {
                    
                }
                
                Matrix() :
                n(0), m(0), size(0), a(nullptr) {
                    
                }
                
                Matrix(T* _a);
                Matrix(const Matrix<T>& M);
                
                template <class T2>
                Matrix(const Matrix<T2>& M);
                
                ~Matrix();
                void reSize(int _n, int _m);
                Matrix<T>& operator=(const Matrix<T>& M);
                
                
                
                //  --- GETS AND SETS ---
                //  Gets
                Matrix<T> getRow(int _i);                                                 //  Devuelve la fila i.
                Matrix<T> getColumn(int _j);                                              //  Devuelve la columna j.
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
                T& operator()(int _k) const;                                            //  Cambia el valor del elemento k. (Notación pseudovectorial).
                
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
                void zero();                                                            //  Llena la matriz de ceros.
                void ones();                                                            //  Llena la matriz de unos.
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
template <class T> inline cda::math::containers::Matrix<T>
operator*(const T& D, const cda::math::containers::Matrix<T>& M);                                    //  Define el producto de un escalar y una matriz por la izquierda.
template <class T> inline cda::math::containers::Matrix<T>
operator&&(const cda::math::containers::Matrix<T>& Mi, const cda::math::containers::Matrix<T>& Md);                           //  Devuelve una matriz que tendrá unidas dos matrices.
template <class T> inline void
operator||(cda::math::containers::Matrix<T>& Mi, cda::math::containers::Matrix<T>& Md);                                       //  Separa una matriz en dos.
template <class T> inline void
operator>>(std::istream& in, cda::math::containers::Matrix<T>& M);                                        //  Importa una matriz desde un fichero.
template <class T> inline std::ostream&
operator<<(std::ostream& out, const cda::math::containers::Matrix<T>& M);                                 //  Exporta la matriz a un ostream.

#endif /* math_containers_matrix_hpp */
