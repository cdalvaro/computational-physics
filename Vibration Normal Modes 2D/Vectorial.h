//
//  Vectorial.h
//  Class for vectorial operations
//  Vibration Normal Modes 2D
//
//  Created by Carlos David on 5/8/12.
//  Copyright (c) 2012 cdalvaro. All rights reserved.
//

#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <string>
#include <stdexcept>

//  Para QR
#define Qmatrix     0x01
#define Rmatrix     0x02
#define QRmatrix    0x04

//  Para mostrar las iteraciones de los autovalores y de los autovectores
#define EVa_Ite     0x01
#define EVe_Ite     0x02
#define ORGANIZED   0x04


namespace CDA {
    
    //  --- VECTOR CLASS ---
    template <class T>
    class Vector {
    private:
        size_t n;
        T *v;
        
    public:
        
        /**
         Class constructor
         
         @param n The size of the vector
         */
        Vector(const size_t &size) :
        n(size), v(nullptr) {
            if (void *mem = std::realloc(v, n * sizeof(T))) {
                v = static_cast<T *>(mem);
            } else {
                throw std::bad_alloc();
            }
        }
        
        /**
         Class default constructor
         */
        Vector() :
        Vector(0) {
            
        }
        
        /**
         Class constructor
         
         @param n The size of the vector
         */
        Vector(const size_t &size, const T &value) :
        n(size), v(nullptr) {
            if (void *mem = std::realloc(v, n * sizeof(T))) {
                v = static_cast<T *>(mem);
            } else {
                throw std::bad_alloc();
            }
            
            std::fill(v, v + n, value);
        }
        
        /**
         Copy constructor
         
         @param vector The source vector
         */
        Vector(const Vector<T> &vector) :
        n(vector.n), v(nullptr) {
            if (void *mem = std::realloc(v, n * sizeof(T))) {
                v = static_cast<T *>(mem);
            } else {
                throw std::bad_alloc();
            }
            
            std::copy(vector.v, vector.v + n, v);
        }
        
        // FIXME: This function is not working properly
        /**
         Move constructor
         
         @param vector The source vector
         */
//        Vector(Vector<T> &&vector) :
//        n(vector.n), v(std::move(vector.v)) {
//
//        }
        
        /**
         Class destructor
         */
        ~Vector() {
            std::free(v);
            v = nullptr;
        }
        
        T *Begin() const {
            return this->v;
        }
        
        T *End() const {
            return this->v + this->n;
        }
        
        /**
         Change the size of the vector

         @param size The new size of the vector
         @param fill If true and size is bigger than the previous size
                     the new extra elements will be set to 0
         */
        void Resize(const size_t &size, const bool &fill = true) {
            if (void *mem = std::realloc(v, size * sizeof(T))) {
                v = static_cast<T *>(mem);
            } else {
                throw std::bad_alloc();
            }
            
            if (fill && size > n) {
                std::fill_n(v + n, size - n, static_cast<T>(0));
            }
            
            n = size;
        }
        
        Vector<T> &operator=(const Vector<T> &vector) {
            if (&vector != this) {
                this->Resize(vector.n, false);
                std::copy(vector.v, vector.v + n, v);
            }
            
            return (* this);
        }
        
        //  FIXME: This function is not working properly
//        Vector<T> &operator=(Vector<T> &&vector) {
//            // TODO: Check method
//            if (&vector != this) {
//                n = vector.n;
//                std::free(v);
//                v = std::move(vector.v);
//            }
//            
//            return (* this);
//        }
        
        /**
         Returns \p elements starting from element \p first_element

         @param first_element The first element to be recovered
         @param elements The number of elements to be recovered
         
         @return A vector with size \p elements and the values starting from \p first_element
         */
        Vector<T> Get(const size_t &first_element, const size_t &elements) {
            if (n - first_element < elements) {
                throw std::out_of_range("There are not enough elements inside the vector");
            }
            
            Vector<T> tmp(n - first_element);
            std::copy(v + first_element, v + first_element + elements, tmp.v);
            
            return tmp;
        }
        
        /**
         Returns a vector starting at position \p first_element

         @param first_element The position of the first element to be recovered
         
         @return A vector with all values starting from \p first_element
         */
        Vector<T> Get(const size_t &first_element) {
            return Get(first_element, n - first_element);
        }
        
        /**
         Copy \p elements from \p vector starting at position \p first_element

         @param first_element The position of the first element to be copied
         @param vector The vector with the elements to be copied
         @param lenght The number of elements to copy
         */
        void Set(const size_t &first_element, const Vector<T> &vector, const size_t &elements) {
            if (first_element + elements > n || elements > vector.n) {
                throw std::out_of_range("Out of bounds");
            }
            
            std::copy(vector.v, vector.v + elements, v + first_element);
        }
        
        /**
         Copy all elements of \p vector at position \p first_element
         
         @param first_element The position of the first element to be copied
         @param vector The vector to copy
         */
        void Set(const size_t &first_element, const Vector<T> &vector) {
            Set(first_element, vector, vector.n);
        }
        
        size_t Size() const {
            return n;
        }
        
        T MaximumElement() const {
            auto max = *v;
            const auto it_end = v + n;
            
            for (auto it = v + 1; it != it_end; ++it) {
                if (*it > max) {
                    max = *it;
                }
            }
            
            return max;
        }
        
        T AbsoluteMaximumElement() const {
            auto max = std::abs(*v);
            const auto it_end = v + n;
            
            T abs_max;
            for (auto it = v + 1; it != it_end; ++it) {
                abs_max = std::abs(*it);
                if (abs_max > max) {
                    max = abs_max;
                }
            }
            
            return max;
        }
        
        T MinimumElement() const {
            auto min = *v;
            const auto it_end = v + n;
            
            for (auto it = v +1; it != it_end; ++it) {
                if (*it < min) {
                    min = *it;
                }
            }
            
            return min;
        }
        
        T AbsoluteMinimumElement() const {
            auto min = std::abs(*v);
            const auto it_end = v + n;
            
            T abs_min;
            for (auto it = v +1; it != it_end; ++it) {
                abs_min = std::abs(*it);
                if (abs_min < min) {
                    min = abs_min;
                }
            }
            
            return min;
        }
        
        T Sum() const {
            auto sum = *v;
            const auto it_end = v + n;
            
            for (auto it = v + 1; it != it_end; ++it) {
                sum += *it;
            }
            
            return sum;
        }
        
        const T &At(const size_t &element) const {
            if (element >= n) {
                throw std::out_of_range("Index out of bounds");
            }
            
            return v[element];
        }
        
        const T &operator[](const size_t &element) const {
            return v[element];
        }
        
        T &At(const size_t &element) {
            if (element >= n) {
                throw std::out_of_range("Index out of bounds");
            }
            
            return v[element];
        }
        
        T &operator[](const size_t &element) {
            return v[element];
        }
        
        T Norm() const {
            return std::sqrt(SquaredNorm());
        }
        
        T SquaredNorm() const {
            return (* this) * (* this);
        }
        
        /**
         Returns a normalized vector

         @return The normalized vector
         */
        Vector<T> Unitary() const {
            return (* this) / Norm();
        }
        
        void Sort() {
            bool change;
            const auto it_end = v + n;
            
            do {
                change = false;
                for (auto it_prev = v, it = it_prev + 1; it != it_end; ++it_prev, ++it) {
                    if (*it_prev > *it) {
                        change = true;
                        std::swap(*it, *it_prev);
                    }
                }
            } while (change);
            
        }
        
        void Fill(const T &value) {
            std::fill(v, v + n, value);
        }
        
        static Vector<T> Zero(const size_t &size) {
            return Vector<T>(size, 0);
        }
        
        void Zero() {
            Fill(0);
        }
        
        static Vector<T> Ones(const size_t &size) {
            return Vector<T>(size, 1);
        }
        
        void Ones() {
            Fill(1);
        }
        
        static Vector<T> Random(const size_t &size) {
            Vector<T> random(size);
            random.Random();
            return random;
        }
        
        void Random(const T &min = 0, const T &max = 1) {
            const auto it_end = v + n;
            for (auto it = v; it != it_end; ++it) {
                *it = static_cast<T>(drand48() * max + min);
            }
        }
        
        bool IsNull() const {
            const auto it_end = v + n;
            for (auto it = v; it != it_end; ++it) {
                if (*it != 0) {
                    return false;
                }
            }
            
            return true;
        }
        
        bool HasDuplicated(const T &precision) const {
            const auto it_end = v + n;
            T distance;
            
            for (auto it1 = v; it1 != it_end; ++it1) {
                for (auto it2 = v; it2 != it_end; ++it2) {
                    if (it1 != it2) {
                        distance = *it1 - *it2;
                        if (std::sqrt(distance * distance) < precision) {
                            return true;
                        }
                    }
                }
            }
            
            return false;
        }
        
        bool HasDuplicated() const {
            const auto it_end = v + n;
            for (auto it1 = v; it1 != it_end; ++it1) {
                for (auto it2 = v; it2 != it_end; ++it2) {
                    if (it1 != it2 && *it1 == *it2) {
                        return true;
                    }
                }
            }
            
            return false;
        }
        
        void Write(const std::string &filename, const std::string &path = getenv("HOME")) const {
            const std::string file_path(path + filename);
            std::ofstream output(file_path.data());
            output << (* this);
        }
        
        //  --- OPERATORS ---
        inline Vector<T> operator+(const Vector<T> &vector) const {
            if (this->n != vector.n) {
                throw std::logic_error("Unable to sum two vector of different size");
            }
            
            Vector<T> new_vector(this->n);
            auto it_new = new_vector.v;
            
            auto it_vector = vector.v;
            const auto it_end = this->v + this->n;
            
            for (auto it_this = this->v; it_this != it_end; ++it_this, ++it_vector, ++it_new) {
                *it_new = *it_this + *it_vector;
            }
            
            return new_vector;
        }
        
        inline Vector<T> &operator+=(const Vector<T> &vector) {
            if (this->n != vector.n) {
                throw std::logic_error("Unable to sum two vector of different size");
            }
            
            auto it_vector = vector.v;
            const auto it_end = this->v + this->n;
            
            for (auto it_this = this->v; it_this != it_end; ++it_this, ++it_vector) {
                *it_this += *it_vector;
            }
            
            return *this;
        }
        
        inline Vector<T> operator-(const Vector<T> &vector) const {
            if (this->n != vector.n) {
                throw std::logic_error("Unable to sum two vector of different sizes");
            }
            
            Vector<T> new_vector(this->n);
            auto it_new = new_vector.v;
            
            auto it_vector = vector.v;
            const auto it_end = this->v + this->n;
            
            for (auto it_this = this->v; it_this != it_end; ++it_this, ++it_vector, ++it_new) {
                *it_new = *it_this - *it_vector;
            }
            
            return new_vector;
        }
        
        inline Vector<T> &operator-=(const Vector<T> &vector) {
            if (this->n != vector.n) {
                throw std::logic_error("Unable to sum two vector of different size");
            }
            
            auto it_vector = vector.v;
            const auto it_end = this->v + this->n;
            
            for (auto it_this = this->v; it_this != it_end; ++it_this, ++it_vector) {
                *it_this -= *it_vector;
            }
            
            return *this;
        }
        
        inline Vector<T> operator*(const T &value) const {
            Vector<T> new_vector(this->n);
            const auto it_end = new_vector.v + new_vector.n;
            
            for (auto it_new = new_vector.v; it_new != it_end; ++it_new) {
                *it_new *= value;
            }
            
            return new_vector;
        }
        
        inline Vector<T> &operator*=(const T &value) {
            const auto it_end = this->v + this->n;
            for (auto it = this->v; it != it_end; ++it) {
                *it *= value;
            }
            
            return *this;
        }
        
        template<class Matrix>
        inline Vector<T> operator*(const Matrix& matrix) const {
            const auto &rows = matrix.Rows();
            if (rows != this->n) {
                throw std::logic_error("This vector and this matrix are incompatible");
            }
            
            const auto &columns = matrix.Columns();
            Vector<T> new_vector(columns, 0);
            
            // TODO: Check this operation
            auto it_new = new_vector.v;
            for (size_t column = 0; column < columns; ++column, ++it_new) {
                auto it_this = this->v;
                for (size_t row = 0; row < rows; ++row, ++it_this) {
                    *it_new += *it_this * matrix(row, column);
                }
            }
            
            return new_vector;
        }
        
        inline Vector<T> operator/(const T &value) const {
            Vector<T> new_vector(this->n);
            auto it_new = new_vector.v;
            const auto it_end = this->v + this->n;
            
            for (auto it = this->v; it != it_end; ++it, ++it_new) {
                *it_new = *it / value;
            }
            
            return new_vector;
        }
        
        inline Vector<T> &operator/=(const T &value) {
            auto it_end = this->v + this->n;
            for (auto it = this->v; it != it_end; ++it) {
                *it /= value;
            }
            
            return *this;
        }
        
        template<typename Integer,
                 typename = std::enable_if<std::is_integral<Integer>::value>>
        inline Vector<T> operator%(const Integer &value) const {
            Vector<Integer> new_vector(this->n);
            auto it_new = new_vector->v;
            const auto it_end = this->v + this->n;
            
            for (auto it = this->v; it != it_end; ++it, ++it_new) {
                *it_new = static_cast<Integer>(*it) % value;
            }
            
            return new_vector;
        }
        
        template<typename Integer,
                 typename = std::enable_if<std::is_integral<Integer>::value>>
        inline Vector<T> &operator%=(const Integer& value) {
            const auto it_end = this->v + this->n;
            for (auto it = this->v; it != it_end; ++it) {
                *it = static_cast<Integer>(*it) % value;
            }
            
            return *this;
        }
        
        inline Vector<T> operator-() const {
            Vector<T> new_vector(this->n);
            auto it_new = new_vector.v;
            const auto it_end = this->v + this->n;
            
            for (auto it = this->v; it != it_end; ++it, ++it_new) {
                *it_new = -(*it);
            }
            
            return new_vector;
        }
        
        Vector<T> CrossProduct3D(const Vector<T> &vector) const {
            if (this->n != vector.n || this->n != 3) {
                throw std::logic_error("Both vectors must be of the same size, and size must be 3");
            }
            
            Vector<T> new_vector(this->n);
            new_vector.v[0] = this->v[1] * vector.v[2] - this->v[2] * vector.v[1];
            new_vector.v[1] = this->v[2] * vector.v[0] - this->v[0] * vector.v[2];
            new_vector.v[2] = this->v[0] * vector.v[1] - this->v[1] * vector.v[0];
            
            return new_vector;
        }
        
        inline T operator*(const Vector<T> &vector) const {
            if (this->n != vector.n) {
                throw std::logic_error("Both vectors must be of the same size");
            }
            
            T value(0);
            auto it_vector = vector.v;
            const auto it_end = this->v + this->n;
            
            for (auto it_this = this->v; it_this != it_end; ++it_this, ++it_vector) {
                value += *it_this * *it_vector;
            }
            
            return value;
        }
        
    };
    
    
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
}

//  --- MORE OPERATORS ---
template <typename T>
inline static CDA::Vector<T> operator*(const T &value, const CDA::Vector<T> &vector) {
    CDA::Vector<T> new_vector(vector.v);
    auto it_new = new_vector.v;
    const auto it_end = vector.v + vector.n;
    
    for (auto it = vector.v; it != it_end; ++it, ++it_new) {
        *it = *it * value;
    }
    
    return new_vector;
}

//  --- MORE OPERATORS ---
template <class T> inline CDA::Matrix<T>
operator*(const T& D, const CDA::Matrix<T>& M);                                    //  Define el producto de un escalar y una matriz por la izquierda.
template <class T> inline CDA::Matrix<T>
operator&&(const CDA::Matrix<T>& Mi, const CDA::Matrix<T>& Md);                           //  Devuelve una matriz que tendrá unidas dos matrices.
template <class T> inline void
operator||(CDA::Matrix<T>& Mi, CDA::Matrix<T>& Md);                                       //  Separa una matriz en dos.
template <class T> inline void
operator>>(std::istream& in, CDA::Matrix<T>& M);                                        //  Importa una matriz desde un fichero.
template <class T> inline std::ostream&
operator<<(std::ostream& out, const CDA::Matrix<T>& M);                                 //  Exporta la matriz a un ostream.
