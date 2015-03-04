//
//  Solve.h
//  Vibration Normal Modes 2D
//
//  Created by Carlos David on 5/8/12.
//  Copyright (c) 2012 NEPTUNO. All rights reserved.
//

#pragma once

//  BOUNDARY CONDITIONS
//  Condiciones en la función
#define BCL_f       0x01    //  BCL - Boundary condition left
#define BCR_f       0x02    //  BCR - Boundary condition right
#define BCT_f       0x04    //  BCT - Boundary condition top
#define BCB_f       0x08    //  BCB - Boundary condition bottom
//  Condiciones en la derivada
#define BCL_df      0x10
#define BCR_df      0x20
#define BCT_df      0x40
#define BCB_df      0x80
//  Condición de contorno interior especial
#define BCI_f       0x10   //  BCI - Boundary condition intern

//  OPTIONS
#define ITERATIONS      0x01    //  Mostrará el número de iteraciones realizadas por el método.
#define SAVE_DATA       0x02
#define IMPORT_DATA     0x04
#define PLOT_MATLAB     0x08    //  Devolverá dos archivos en la ruta indicada, un .csv con los datos
                                //  y un .m con las instrucciones para pintar en MATLAB
#define DESKTOP         0x10
#define DOCUMENTS       0x20

//  Para el método de integración de diferencias finitas
#define LUmethod        0x20
#define GSmethod        0x40


#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#import "Vectorial.cpp"

typedef double EDP_T;
#define vectorEDP_T vector<EDP_T>
#define matrixEDP_T matrix<EDP_T>

using namespace std;

class EDP {
private:
    //  ECUACIÓN DE ONDAS
    //  Para almacenar la situación anterior.
    vectorEDP_T old1D;
    matrixEDP_T old2D;
    matrix<bool> fixedEDP;
    
public:
    //  --- RUTA PARA GUARDAR DATOS ---
    string pathEDP;
    
    
    //  --- ECUACIONES DIFERENCIALES EN DERIVADAS PARCIALES ---
    
    
    //  -- MÉTODO DE LAS DIFERENCIAS FINITAS - SISTEMAS DE 1 DIMENSIÓN --
    //  Resuelve ecuaciones diferenciales del tipo: Y''(x) + A(x)·Y'(x) + B(x)·Y(x) = C(x) + D
    
    //  PARÁMETROS DE DEFINICIÓN DE LA FUNCIÓN
    EDP_T (* A)(EDP_T x), (* B)(EDP_T x), (* C)(EDP_T x);     //  C debe definirse como: C(x) + D
    
    vectorEDP_T solveDIF_FIN(unsigned char bc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, EDP_T err);
    
    
    
    //  CONDICIONES DE CONTORNO NORMALES -> PARA TODAS LAS FUNCIONES
    EDP_T (* BCL)(EDP_T x, EDP_T y), (* BCR)(EDP_T x, EDP_T y), (* BCT)(EDP_T x, EDP_T y), (* BCB)(EDP_T x, EDP_T y);
    
    //  CONDICIONES DE CONTORNO ESPECIALES -> LAPLACE Y POISSON
    EDP_T (* SBCL)(EDP_T uR, EDP_T uT, EDP_T uB);
    EDP_T (* SBCR)(EDP_T uL, EDP_T uT, EDP_T uB);
    EDP_T (* SBCT)(EDP_T uL, EDP_T uR, EDP_T uB);
    EDP_T (* SBCB)(EDP_T uL, EDP_T uR, EDP_T uT);
    EDP_T (* SBCI)(EDP_T uL, EDP_T uR, EDP_T uT, EDP_T uB);
    
    
    //  -- ECUACIONES DE LAPLACE Y POISSON - SISTEMAS DE 2 DIMENSIONES --
    //  Resuelve las ecuaciones diferenciales: ∆U(x,y) = 0 y ∆U(x,y) = F, F ≠ F(U)
    
    //  - ECUACIÓN DE LAPLACE
    //  Función principal / Método de resolución
    matrixEDP_T solveLAPLACE(unsigned char bc, unsigned char sbc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err, int tol);
    
    //  Funciones de llamada.
    matrixEDP_T solveLAPLACE(unsigned char bc, unsigned char sbc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err);
    matrixEDP_T solveLAPLACE(unsigned char bc, unsigned char sbc, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err, int tol);
    matrixEDP_T solveLAPLACE_NBC(unsigned char bc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err);
    matrixEDP_T solveLAPLACE_SBC(unsigned char sbc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err);
    matrixEDP_T solveLAPLACE_NBC(unsigned char bc, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err);
    matrixEDP_T solveLAPLACE_SBC(unsigned char sbc, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err);
    
    
    //  POISSON
    //  CONDICIÓN DE LA EC. DE POISSON
    //  Es importante que esta función no dependa de otros puntos de la matriz.
    EDP_T (* F)(EDP_T x, EDP_T y);
    
    //  - ECUACIÓN DE POISSON
    //  Función principal / Método de resolución
    matrixEDP_T solvePOISSON(unsigned char bc, unsigned char sbc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err, int tol);
    
    //  Funciones de llamada
    matrixEDP_T solvePOISSON(unsigned char bc, unsigned char sbc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err);
    matrixEDP_T solvePOISSON(unsigned char bc, unsigned char sbc, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err, int tol);
    matrixEDP_T solvePOISSON_NBC(unsigned char bc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err);
    matrixEDP_T solvePOISSON_SBC(unsigned char sbc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err);
    matrixEDP_T solvePOISSON_NBC(unsigned char bc, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err);
    matrixEDP_T solvePOISSON_SBC(unsigned char sbc, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err);
    
    
    
    //  -- ECUACIONES DE ONDA Y DE CALOR --
    //  Resuelve las EDPs de onda: ∂²u/∂t² = Q()·∆u y de calor: ∂u/∂t = Q()·∆u
    
    //  PARÁMETROS
    EDP_T time, dt;
    EDP_T (* Q1D)(EDP_T x);  //  Constante o función Q() que acompaña al laplaciano en 1 dimensión.
    EDP_T (* Q2D)(EDP_T x, EDP_T y);
    
    //  FUNCIONES DE LA EC. DE ONDA
    vectorEDP_T solveWAVE(unsigned char bc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& cI, vectorEDP_T& cId);
    vectorEDP_T solveWAVE(unsigned char bc, vectorEDP_T& x, vectorEDP_T& cI, vectorEDP_T& cId);
    matrixEDP_T solveWAVE(unsigned char bc, unsigned char opt, vectorEDP_T &x, vectorEDP_T &y, matrixEDP_T &cI, matrixEDP_T &cId, matrix<bool> &fixed);
    matrixEDP_T solveWave(unsigned char bc, vectorEDP_T &x, vectorEDP_T &y, matrixEDP_T &cI, matrixEDP_T &cId, matrix<bool> &fixed);
    matrixEDP_T solveWAVE(unsigned char bc, vectorEDP_T &x, vectorEDP_T &y, matrixEDP_T &cI, matrixEDP_T &cId);
    
    
    //  FUNCIONES DE LA EC. DE CALOR
    //  Elementos necesarios para la función
    EDP_T theta;
    vectorEDP_T solveHEAT(unsigned char bc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y);
    vectorEDP_T solveHEAT(unsigned char bc, vectorEDP_T& x, vectorEDP_T& y);
    matrixEDP_T solveHEAT(unsigned char bc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI);
    matrixEDP_T solveHEAT(unsigned char bc, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI);
    
    
    
    //  -- ECUACIONES DEL TIPO -> ∂²u/∂t² = -k²u --
    //  Calcula los autovalores y autovectores de dicha ecuación.
    vectorEDP_T eigenVAL_VEC(vectorEDP_T& x, int mode, unsigned char bc,unsigned char opt);
    vectorEDP_T eigenVAL_VEC(vectorEDP_T& x, int mode, unsigned char bc);
    matrixEDP_T eigenVAL_VEC(vectorEDP_T &x, vectorEDP_T &y, int modeX, int modeY, unsigned char bc, unsigned char opt);
    matrixEDP_T eigenVAL_VEC(vectorEDP_T &x, vectorEDP_T &y, int modeX, int modeY, unsigned char bc);
    
    
    //  -- FUNCIONES DE SALIDA --
    void saveDATA(const string path, const string fileName, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& forSave);
    void saveDATA(const string path, const string fileName, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& forSave);
    void saveDATA(const string fileName, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& forSave);
    
    void printMATLAB(const string path, const string fileName, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& forPrint);
    void printMATLAB(const string fileName, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& forPrint);
};
