//
//  Vectorial.h
//  Class for vectorial operations
//  Vibration Normal Modes 2D
//
//  Created by Carlos David on 5/8/12.
//  Copyright (c) 2012 NEPTUNO. All rights reserved.
//

#pragma once

#define MAT_TEMPLATE template <class T>
#define vectorT vector<T>
#define matrixT matrix<T>

//  Para QR
#define Qmatrix     0x01
#define Rmatrix     0x02
#define QRmatrix    0x04

//  Para mostrar las iteraciones de los autovalores y de los autovectores
#define EVa_Ite     0x01
#define EVe_Ite     0x02
#define ORGANIZED   0x04

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <string>

using namespace std;

MAT_TEMPLATE
class vector;

MAT_TEMPLATE
class matrix;


                //  --- VECTOR CLASS ---
MAT_TEMPLATE
class vector {
private:
    int n;
    T* v;
    
public:
    //  --- DEFINITION ---
    vector(int _n);
    vector();
    vector(T* _v);
    vector(const vectorT& V);
    ~vector();
    void reSize(int _n);
    vectorT& operator=(const vectorT& V);
    
    
    
    //  --- GETS AND SETS ---
    //  Gets
    vectorT get(int _i, int lenght);                                        //  Devuelve un vector de longitud lenght con los elementos desde el elemento i.
    vectorT get(int _i);                                                    //  Devuelve un vector con los elementos desde el elemento i.
    
    //  Sets
    void set(int _i, const vectorT& V, int lenght);                         //  Permite introducir los elementos de un vector en otro.
    void set(int _i, const vectorT& V);                                     //  Permite introducir los elementos de un vector en otro.
    
    
    
    //  --- FUNCTIONS ---
    //  Resturns an integer / a float / a double
    int dim() const;                                                        //  Devuelve la dimensión del vector.
    T max() const;                                                          //  Devuelve el mayor elemento del vector.
    T min() const;                                                          //  Devuelve el menor elemento del vector.
    T maxAbs() const;                                                       //  Devuelve el mayor elemento del vector en valor absoluto.
    T sum() const;                                                          //  Suma todos los elementos del vector.
    T& operator()(int _i) const;                                            //  Devuelve el elemento i del vector
    
    //  Returns a vector / a matrix
    vectorT mod(const T& D);                                                //  Genera un vector de módulo D.
    vectorT unitary();                                                      //  Normaliza el vector.
    vectorT organize();                                                     //  Ordena los elementos del vector del mínimo al máximo.
    matrixT transpose();                                                    //  Transpone un vector devolviendo una matriz columna.
    vectorT zero(int _n);                                                   //  Genera un vector de ceros de dimensión n.
    vectorT ones(int _n);                                                   //  Genera un vector de unos de dimensión n.
    
    //  Returns a bool
    bool null();                                                            //  Comprueba si el vector es un vector de ceros.
    bool duplicate(T range);                                                //  Comprueba si hay elementos duplicados (próximos);
    bool duplicate();                                                       //  Comprueba si hay elementos duplicados.
    
    //  Void functions
    void ones();                                                            //  Llena el vector de unos.
    void zero();                                                            //  LLena el vector de ceros.
    void rand();                                                            //  Llena el vector de números aleatorios.
    void rand(const T& Min, const T& Max);                                  //  LLena el vector de números aleatorios decimales entre dos valores.
    void write(const string& path, const string& filename);                 //  Escribe el vector en un fichero.
    void write(const string& filename);
    
    
    
    //  --- OPERATORS ---
    vectorT operator+(const vectorT& V);                                    //  Definición de la suma de dos vectores.
    vectorT& operator+=(const vectorT& V);                                  //  Suma del vector con otro sobre sí mismo
    vectorT operator-(const vectorT& V);                                    //  Definición de la diferencia de dos vectores.
    vectorT& operator-=(const vectorT& V);                                  //  Diferencia del vector con otro sobre sí mismo.
    vectorT operator*(const T& D);                                          //  Producto del vector por un escalar por la derecha.
    vectorT& operator*=(const T& D);                                        //  Producto del vector por un esaclar sobre sí mismo.
    vectorT operator*(const matrixT& Mc);                                   //  Producto del vector por una matriz.
    vectorT operator/(const T& D);                                          //  Cociente del vector por un escalar.
    vectorT& operator/=(const T& D);                                        //  Cociente del vector por un escalar sobre sí mismo.
    vectorT operator%(const int& I);                                        //  Resto de la división del vector por un número.
    vectorT& operator%=(const int& I);                                      //  Resto de la división del vector por un número sobre sí mismo.
    vectorT operator-();                                                    //  Cambio de signo de los elementos del vector.
    vectorT cross(const vectorT& V);                                        //  Definición del producto vectorial de dos vectores.
    T operator!();                                                          //  Norma del vector.
    T operator~();                                                          //  Norma cuadrática del vector.
    T operator*(const vectorT& V);                                          //  Definición del producto escalar.
};

//  --- MORE OPERATORS ---
MAT_TEMPLATE inline
vectorT operator*(const T& D, const vectorT& V);                            //  Exporta el vector a un ostream.
MAT_TEMPLATE inline void
operator>>(istream& in, vectorT& V);                                        //  Importa un vector desde un fichero.
MAT_TEMPLATE inline ostream&
operator<<(ostream& out, const vectorT& V);                                 //  Define el producto de un escalar y un vector por la izquierda.



//  --- OTHER FUNCTIONS ---
MAT_TEMPLATE inline vectorT
sqrt(const vectorT& V);                                                     //  Calcula la raíz cuadrada de los elementos del vector.
MAT_TEMPLATE inline vectorT
pow(const vectorT& V, double exp);                                          //  Eleva todos los elementos del vector a una potencia.
MAT_TEMPLATE inline vectorT
log(const vectorT& V);                                                      //  Calcula el logaritmo de cada uno de los elementos del vector.
MAT_TEMPLATE inline vectorT
zero(int _n);                                                               //  Genera un vector de ceros de dimensión n.
MAT_TEMPLATE inline vectorT
ones(int _n);                                                               //  Genera un vector de unos de dimensión n.
MAT_TEMPLATE inline vectorT
rand(int _n);                                                               //  Llena el vector de números aleatorios.
MAT_TEMPLATE inline vectorT
rand(int _n, const T& Min, const T& Max);                                   //  LLena el vector de números aleatorios decimales entre dos valores.
MAT_TEMPLATE inline vectorT
round(const vectorT& V);                                                    //  Redondea los valores de un vector.



                //  --- MATRIX CLASS ---
MAT_TEMPLATE
class matrix {
private:
    int n, m, size;
    T* a;
    
public:
    //  --- DEFINITION ---
    matrix(int _n, int _m);
    matrix();
    matrix(T* _a);
    matrix(const matrixT& M);
    
    template <class T2>
    matrix(const matrix<T2>& M);
    
    ~matrix();
    void reSize(int _n, int _m);
    matrixT& operator=(const matrixT& M);
    
    
    
    //  --- GETS AND SETS ---
    //  Gets
    matrixT getRow(int _i);                                                 //  Devuelve la fila i.
    matrixT getColumn(int _j);                                              //  Devuelve la columna j.
    matrixT getMatrix(int _i, int _j, int rows, int columns);               //  Devuelve la submatriz (i,j) con tamaño (n,m).
    matrixT getMatrix(int _i, int _j);                                      //  Devuelve la submatriz (i,j).
    
    vectorT getDiagonal();                                                  //  Devuelve la diagonal de una matriz.
    vectorT getRowV(int _i, int _j);                                        //  Devuelve como vector la fila i desde la columna j.
    vectorT getRowV(int _i);                                                //  Devuelve como vector la fila i.
    vectorT getColumnV(int _i, int _j);                                     //  Devuelve como vector la columna j desde la fila i.
    vectorT getColumnV(int _j);                                             //  Devuelve como vector la columna j.
    
    //  Sets
    void setColumn(int _i, int _j, const matrixT& Mc);                      //  Introduce una matriz columna en la columna j desde la fila i.
    void setColumn(int _j, const matrixT& Mc);                              //  Introduce una matriz columna en la columna j.
    void setColumnV(int _i, int _j, const vectorT& V);                      //  Introduce un vector en la columna j desde la fila i.
    void setColumnV(int _j, const vectorT& V);                              //  Introduce un vector en la columna j.
    void setRow(int _i, int _j, const matrixT& Mr);                         //  Introduce una matriz fila en la fila i desde la columna j.
    void setRow(int _i, const matrixT& Mr);                                 //  Introduce una matriz fila en la fila i.
    void setRowV(int _i, int _j, const vectorT& V);                         //  Introduce un vector en la fila i desde la columna j.
    void setRowV(int _i, const vectorT& V);                                 //  Introduce un vector en la fila i.
    void setMatrix(int _i, int _j, const matrixT& M);                       //  Introduce una matriz en la posición (i,j).
    
    
    
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
    
    //  Returns a matrix / a vector
    matrixT sumRows();                                                      //  Suma todos los elementos de las filas y las devuelve en una matriz columna.
    matrixT sumColumns();                                                   //  Suma todos los elementos de las columnas y las devuelve en una matriz fila.
    matrixT transpose();                                                    //  Transpone la matriz.
    matrixT zero(int _n, int _m);                                           //  Devuelve una matriz de tamaño (n,m) llena de ceros.
    matrixT ones(int _n, int _m);                                           //  Devuelve una matriz de tamaño (n,m) llena de unos.
    matrixT identity(int _n, int _m);                                       //  Devuelve una matriz identidad de tamaño (n,m).
    
    vectorT sumRowsV();                                                     //  Suma todos los elementos de las filas y las devuelve en un vector.
    vectorT sumColumnsV();                                                  //  Suma todos los elementos de las columnas y las devuelve en un vector.
    
    //  Returns a bool
    bool null();                                                            //  Comprueba si la matriz es una matriz de ceros.
    bool duplicate(T range);                                                //  Comprueba si hay elementos duplicados (próximos);
    bool duplicate();                                                       //  Comprueba si hay elementos duplicados.
    
    //  Void functions
    void zero();                                                            //  Llena la matriz de ceros.
    void ones();                                                            //  Llena la matriz de unos.
    void identity();                                                        //  Combierte a la matriz en la identidad.
    void write(const string& path, const string& filename);                 //  Escribe la matriz en un fichero.
    void write(const string& filename);
    
    
    
    //  --- MÉTODO LU ---
    matrixT LU();                                                           //  Descompone una matriz en dos por el método LU.
    matrixT U();                                                            //  Devuelve la matriz U del método LU.
    matrixT L();                                                            //  Devuelve la matriz L del método LU.
    vectorT solveLU(const vectorT& B);                                      //  Resuelve un sistema de ecuaciones por el método LU.
    vectorT solveLU3d(const vectorT& B);                                    //  Resuelve un sistema de ecuaciones tridiagonal por el método LU.
    vectorT solveGS3d(const vectorT& B, T err);                             //  Resuelve un sistema de ecuaciones tridiagonal por el método Gauss-Seidel.
    
    //  --- MÉTODO QR ---
    matrixT QR(unsigned char QorR);                                         //  Descompone una matriz en dos por el método QR.
    matrixT eigenVectors(int maxIte, unsigned char opt);                    //  Calcula los autovectores a partir de los autovalores (QR).
    matrixT eigenVectors(unsigned char opt);                                //  Calcula los autovectores a partir de los autovaleres (QR).
    matrixT eigenVectors();                                                 //  Calcula los autovectores a partir de los autovalores (QR).
    vectorT eigenVector(T eigenValue, int maxIte, unsigned char opt);       //  Calcula el autovector correspondiente al autovector introducido.
    vectorT eigenVector(T eigenValue, unsigned char opt);                   //  Calcula el autovector correspondiente al autovector introducido.
    vectorT eigenVector(T eigenValue);                                      //  Calcula el autovector correspondiente al autovector introducido.
    vectorT eigenValues(int maxIte, T factor, unsigned char opt);           //  Calcula los autovalores de una matriz por el método QR.
    vectorT eigenValues(int maxIte, unsigned char opt);                     //  Calcula los autovalores de una matriz por el método QR.
    vectorT eigenValues(unsigned char opt);                                 //  Calcula los autovalores de una matriz por el método QR.
    vectorT eigenValues();                                                  //  Calcula los autovalores de una matraz por el método QR.
    
    
    
    //  --- OPERATORS ---
    matrixT& operator=(const T* Array);                                     //  Introduce los elementos de una matriz a través de un array.
    matrixT operator+(const matrixT& M);                                    //  Definición de suma de dos matrices.
    matrixT& operator+=(const matrixT& M);                                  //  Suma la matriz con otra sobre sí misma.
    matrixT operator-(const matrixT& M);                                    //  Definición de diferencia de dos matrices.
    matrixT& operator-=(const matrixT& M);                                  //  Diferencia de la matriz con otra sobre sí misma.
    matrixT operator*(const T& D);                                          //  Producto de la matriz por una constante.
    matrixT& operator*=(const T& D);                                        //  Producto de la matriz sobre sí misma por una constante.
    matrixT operator/(const T& D);                                          //  Cociente de la matriz por una constante.
    matrixT& operator/=(const T& D);                                        //  Cociente de la matriz sobre sí misma por una constante.
    matrixT operator-();                                                    //  Cambia el signo de todos los elementos de la matriz.
    matrixT operator*(const matrixT& M);                                    //  Definición del producto de dos matrices.
    matrixT& operator*=(const matrixT& M);                                  //  Producto de una matriz por otra sobre sí misma.
    matrixT operator*(const vectorT& V);                                    //  Producto de una matriz por un vector.
    matrixT operator^(const int exp);                                       //  Eleva a una potencia la matriz ó calcula su inversa.
};

//  --- MORE OPERATORS ---
MAT_TEMPLATE inline matrixT
operator*(const T& D, const matrixT& M);                                    //  Define el producto de un escalar y una matriz por la izquierda.
MAT_TEMPLATE inline matrixT
operator&&(const matrixT& Mi, const matrixT& Md);                           //  Devuelve una matriz que tendrá unidas dos matrices.
MAT_TEMPLATE inline void
operator||(matrixT& Mi, matrixT& Md);                                       //  Separa una matriz en dos.
MAT_TEMPLATE inline void
operator>>(istream& in, matrixT& M);                                        //  Importa una matriz desde un fichero.
MAT_TEMPLATE inline ostream&
operator<<(ostream& out, const matrixT& M);                                 //  Exporta la matriz a un ostream.



//  --- OTHER FUNCTIONS ---
MAT_TEMPLATE inline matrixT
zero(int _n, int _m);                                                       //  Devuelve una matriz de tamaño (n,m) llena de ceros.
MAT_TEMPLATE inline matrixT
ones(int _n, int _m);                                                       //  Devuelve una matriz de tamaño (n,m) llena de unos.
MAT_TEMPLATE inline matrixT
identity(int _n, int _m);                                                   //  Devuelve una matriz identidad de tamaño (n,m).
MAT_TEMPLATE inline matrixT
setdiff(matrixT& A, const matrixT& B, const int& reps);                     //  Devuelve los valores de A que no se encuentran en B.
MAT_TEMPLATE inline matrixT
setdiff(matrixT& A, const matrixT& B);                                      //  Devuelve los valores de A que no se encuentran en B.
