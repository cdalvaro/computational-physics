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

#include "../containers/vector.hpp"
#include "../containers/matrix.hpp"

namespace cda {
    namespace math {
        namespace differential_equations {

            typedef double EDP_T;
            
            class EDP {
            private:
                
                //  ECUACIÓN DE ONDAS
                //  Para almacenar la situación anterior.
                containers::Vector<EDP_T> old1D;
                containers::Matrix<EDP_T> old2D;
                containers::Matrix<bool> fixedEDP;
                
            public:
                //  --- RUTA PARA GUARDAR DATOS ---
                std::string pathEDP;
                
                //  --- ECUACIONES DIFERENCIALES EN DERIVADAS PARCIALES ---
                
                //  -- MÉTODO DE LAS DIFERENCIAS FINITAS - SISTEMAS DE 1 DIMENSIÓN --
                //  Resuelve ecuaciones diferenciales del tipo: Y''(x) + A(x)·Y'(x) + B(x)·Y(x) = C(x) + D
                
                //  PARÁMETROS DE DEFINICIÓN DE LA FUNCIÓN
                EDP_T (* A)(EDP_T x), (* B)(EDP_T x), (* C)(EDP_T x);     //  C debe definirse como: C(x) + D
                
                containers::Vector<EDP_T> solveDIF_FIN(unsigned char bc, unsigned char opt,
                                                       containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y, EDP_T err);
                
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
                containers::Matrix<EDP_T> solveLAPLACE(unsigned char bc, unsigned char sbc, unsigned char opt,
                                                       containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y,
                                                       containers::Matrix<EDP_T>& cI, EDP_T err, int tol);
                
                //  Funciones de llamada.
                containers::Matrix<EDP_T> solveLAPLACE(unsigned char bc, unsigned char sbc, unsigned char opt,
                                                       containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y,
                                                       containers::Matrix<EDP_T>& cI, EDP_T err);
                
                containers::Matrix<EDP_T> solveLAPLACE(unsigned char bc, unsigned char sbc,
                                                       containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y,
                                                       containers::Matrix<EDP_T>& cI, EDP_T err, int tol);
                
                containers::Matrix<EDP_T> solveLAPLACE_NBC(unsigned char bc, unsigned char opt,
                                                           containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y,
                                                           containers::Matrix<EDP_T>& cI, EDP_T err);
                
                containers::Matrix<EDP_T> solveLAPLACE_SBC(unsigned char sbc, unsigned char opt,
                                                           containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y,
                                                           containers::Matrix<EDP_T>& cI, EDP_T err);
                
                containers::Matrix<EDP_T> solveLAPLACE_NBC(unsigned char bc,
                                                           containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y,
                                                           containers::Matrix<EDP_T>& cI, EDP_T err);
                
                containers::Matrix<EDP_T> solveLAPLACE_SBC(unsigned char sbc,
                                                           containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y,
                                                           containers::Matrix<EDP_T>& cI, EDP_T err);
                
                //  POISSON
                //  CONDICIÓN DE LA EC. DE POISSON
                //  Es importante que esta función no dependa de otros puntos de la matriz.
                EDP_T (* F)(EDP_T x, EDP_T y);
                
                //  - ECUACIÓN DE POISSON
                //  Función principal / Método de resolución
                containers::Matrix<EDP_T> solvePOISSON(unsigned char bc, unsigned char sbc, unsigned char opt,
                                                       containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y,
                                                       containers::Matrix<EDP_T>& cI, EDP_T err, int tol);
                
                //  Funciones de llamada
                containers::Matrix<EDP_T> solvePOISSON(unsigned char bc, unsigned char sbc, unsigned char opt,
                                                       containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y,
                                                       containers::Matrix<EDP_T>& cI, EDP_T err);
                
                containers::Matrix<EDP_T> solvePOISSON(unsigned char bc, unsigned char sbc,
                                                       containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y,
                                                       containers::Matrix<EDP_T>& cI, EDP_T err, int tol);
                
                containers::Matrix<EDP_T> solvePOISSON_NBC(unsigned char bc, unsigned char opt,
                                                           containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y,
                                                           containers::Matrix<EDP_T>& cI, EDP_T err);
                
                containers::Matrix<EDP_T> solvePOISSON_SBC(unsigned char sbc, unsigned char opt,
                                                           containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y,
                                                           containers::Matrix<EDP_T>& cI, EDP_T err);
                
                containers::Matrix<EDP_T> solvePOISSON_NBC(unsigned char bc,
                                                           containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y,
                                                           containers::Matrix<EDP_T>& cI, EDP_T err);
                
                containers::Matrix<EDP_T> solvePOISSON_SBC(unsigned char sbc,
                                                           containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y,
                                                           containers::Matrix<EDP_T>& cI, EDP_T err);
                
                //  -- ECUACIONES DE ONDA Y DE CALOR --
                //  Resuelve las EDPs de onda: ∂²u/∂t² = Q()·∆u y de calor: ∂u/∂t = Q()·∆u
                
                //  PARÁMETROS
                EDP_T time, dt;
                EDP_T (* Q1D)(EDP_T x);  //  Constante o función Q() que acompaña al laplaciano en 1 dimensión.
                EDP_T (* Q2D)(EDP_T x, EDP_T y);
                
                //  FUNCIONES DE LA EC. DE ONDA
                containers::Vector<EDP_T> solveWAVE(unsigned char bc, unsigned char opt,
                                                    containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& cI,
                                                    containers::Vector<EDP_T>& cId);
                
                containers::Vector<EDP_T> solveWAVE(unsigned char bc,
                                                    containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& cI,
                                                    containers::Vector<EDP_T>& cId);
                
                containers::Matrix<EDP_T> solveWAVE(unsigned char bc, unsigned char opt,
                                                    containers::Vector<EDP_T> &x, containers::Vector<EDP_T> &y,
                                                    containers::Matrix<EDP_T> &cI, containers::Matrix<EDP_T> &cId,
                                                    containers::Matrix<bool> &fixed);
                
                containers::Matrix<EDP_T> solveWave(unsigned char bc,
                                                    containers::Vector<EDP_T> &x, containers::Vector<EDP_T> &y,
                                                    containers::Matrix<EDP_T> &cI, containers::Matrix<EDP_T> &cId,
                                                    containers::Matrix<bool> &fixed);
                
                containers::Matrix<EDP_T> solveWAVE(unsigned char bc,
                                                    containers::Vector<EDP_T> &x, containers::Vector<EDP_T> &y,
                                                    containers::Matrix<EDP_T> &cI, containers::Matrix<EDP_T> &cId);
                
                //  FUNCIONES DE LA EC. DE CALOR
                //  Elementos necesarios para la función
                EDP_T theta;
                containers::Vector<EDP_T> solveHEAT(unsigned char bc, unsigned char opt,
                                                    containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y);
                containers::Vector<EDP_T> solveHEAT(unsigned char bc,
                                                    containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y);
                containers::Matrix<EDP_T> solveHEAT(unsigned char bc, unsigned char opt,
                                                    containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y, containers::Matrix<EDP_T>& cI);
                containers::Matrix<EDP_T> solveHEAT(unsigned char bc,
                                                    containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y, containers::Matrix<EDP_T>& cI);
                
                //  -- ECUACIONES DEL TIPO -> ∂²u/∂t² = -k²u --
                //  Calcula los autovalores y autovectores de dicha ecuación.
                containers::Vector<EDP_T> eigenVAL_VEC(containers::Vector<EDP_T>& x,
                                                       int mode, unsigned char bc,unsigned char opt);
                containers::Vector<EDP_T> eigenVAL_VEC(containers::Vector<EDP_T>& x,
                                                       int mode, unsigned char bc);
                containers::Matrix<EDP_T> eigenVAL_VEC(containers::Vector<EDP_T> &x, containers::Vector<EDP_T> &y,
                                                       int modeX, int modeY, unsigned char bc, unsigned char opt);
                containers::Matrix<EDP_T> eigenVAL_VEC(containers::Vector<EDP_T> &x, containers::Vector<EDP_T> &y,
                                                       int modeX, int modeY, unsigned char bc);
                
                //  -- FUNCIONES DE SALIDA --
                void saveDATA(const std::string path, const std::string fileName, unsigned char opt,
                              containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y, containers::Matrix<EDP_T>& forSave);
                void saveDATA(const std::string path, const std::string fileName,
                              containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y, containers::Matrix<EDP_T>& forSave);
                void saveDATA(const std::string fileName, unsigned char opt,
                              containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y,containers:: Matrix<EDP_T>& forSave);
                
                void printMATLAB(const std::string path, const std::string fileName, unsigned char opt,
                                 containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y, containers::Matrix<EDP_T>& forPrint);
                void printMATLAB(const std::string fileName, unsigned char opt,
                                 containers::Vector<EDP_T>& x, containers::Vector<EDP_T>& y, containers::Matrix<EDP_T>& forPrint);
            };
        
        } /* namespace math */
    } /* namespace differential_equations */
} /* namespace cda */
