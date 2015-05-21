//
//  Solve.cpp
//  Vibration Normal Modes 2D
//
//  Created by Carlos David on 5/8/12.
//  Copyright (c) 2012 NEPTUNO. All rights reserved.
//

#include "SolveEDP.h"


int tolDef = 1E+05;     //  Número máximo de iteraciones por defecto.
int errDef = 1E-04;     //  Variación máxima que puede haber entre dos iteraciones.
bool initEDP = false;
const string EDPwarning = "\n[EDP::";


//  --- ECUACIONES DIFERENCIALES EN DERIVADAS PARCIALES ---

//  -- MÉTODO DE LAS DIFERENCIAS FINITAS - SISTEMAS DE 1 DIMENSIÓN --
//  Resuelve ecuaciones diferenciales del tipo: Y''(x) + A(x)·Y'(x) + B(x)·Y(x) = C(x) + D

vectorEDP_T EDP::solveDIF_FIN(unsigned char bc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, EDP_T err)
{
    int dim = x.dim();
    vectorEDP_T sol(dim);
    EDP_T h = (x(dim-1) - x(0))/(EDP_T)(dim-1);
    EDP_T tI = x(0);
    
    x.reSize(x.dim());
    
    if ((bc& BCL_f) != 0 and (bc& BCR_f) != 0)        //  Ambas condiciones de contorno en la función
    {
        matrixEDP_T diagMA(dim-2,3);
        vectorEDP_T b(dim-2), solV(dim-2);
        
        for (int i=0; i<dim-2; i++) {
            tI += h;
            diagMA(i,0) = -((h*A(tI))/2.0 + 1.0);
            diagMA(i,1) = 2.0 - h*h*B(tI);
            diagMA(i,2) = ((h*A(tI))/2.0) - 1.0;
            b(i) = -h*h*C(tI);
        }
        b(0) += (1.0 + (A(x(0)+h)*h)/2.0)*BCL(x(0),0.0);
        b(dim-3) += (1.0 - (A(x(dim-1)-h)*h)/2.0)*BCR(x(dim-1),0.0);
        
        if ((opt& LUmethod) != 0) {
            solV = diagMA.solveLU3d(b);
        } else if ((opt& GSmethod) != 0) {
            solV = diagMA.solveGS3d(b, err);
        } else {
            solV = diagMA.solveLU3d(b);
        }
        
        sol = y;
        sol.set(1, solV);
    }
    
    if ((opt& BCL_df) == BCL_df and (opt& BCR_df) == BCR_df)      //  Ambas condiciones de contorno en la derivada
    {
        matrixEDP_T diagMA(dim,3);
        vectorEDP_T b(dim), solV(dim);
        
        diagMA(0,0) = 0.0;
        diagMA(0,1) = 2.0 - h*h*B(x(0));
        diagMA(0,2) = -2.0;
        b(0) = -h*h*C(x(0)) + 2.0*h*BCL(x(0),0.0)*(h*A(x(0))/2.0 - 1.0);
        for (int i=1; i<dim-1; i++) {
            tI += h;
            diagMA(i,0) = -((h*A(tI))/2.0 + 1.0);
            diagMA(i,1) = 2.0 - h*h*B(tI);
            diagMA(i,2) = ((h*A(tI))/2.0) - 1.0;
            b(i) = -h*h*C(tI);
        }
        diagMA(dim-1,0) = -2.0;
        diagMA(dim-1,1) = 2.0 - h*h*B(x(dim-1));
        diagMA(dim-1,2) = 0.0;
        b(dim-1) = -h*h*C(x(dim-1)) + 2.0*h*BCR(x(dim-1),0.0)*(h*A(x(dim-1))/2.0 + 1.0);
        
        if ((opt& LUmethod) != 0) {
            solV = diagMA.solveLU3d(b);
        } else if ((opt& GSmethod) != 0) {
            solV = diagMA.solveGS3d(b, err);
        } else {
            solV = diagMA.solveLU3d(b);
        }
        
        sol = solV;
    }
    
    if ((opt& BCL_f) == BCL_f and (opt& BCR_df) == BCR_df)       //  Condición de extremo izquierdo en la función y extremo derecho en la derivada
    {
        matrixEDP_T diagMA(dim-1,3);
        vectorEDP_T b(dim-1), solV(dim-1);
        
        for (int i=0; i<dim-2; i++) {
            tI += h;
            diagMA(i,0) = -((h*A(tI))/2.0 + 1.0);
            diagMA(i,1) = 2.0 - h*h*B(tI);
            diagMA(i,2) = ((h*A(tI))/2.0) - 1.0;
            b(i) = -h*h*C(tI);
        }
        b(0) += (1.0 + (A(x(0)+h)*h)/2.0)*BCL(x(0),0.0);
        diagMA(dim-2,0) = -2.0;
        diagMA(dim-2,1) = 2.0 - h*h*B(x(dim-1));
        diagMA(dim-2,2) = 0.0;
        b(dim-2) = -h*h*C(x(dim-1)) + 2.0*h*BCR(x(dim-1),0.0)*(h*A(x(dim-1))/2.0 + 1.0);
        
        if ((opt& LUmethod) != 0) {
            solV = diagMA.solveLU3d(b);
        } else if ((opt& GSmethod) != 0) {
            solV = diagMA.solveGS3d(b, err);
        } else {
            solV = diagMA.solveLU3d(b);
        }
        
        sol = y;
        sol.set(1, solV);
    }
    
    if ((opt& BCL_df) == BCL_df and (opt& BCR_f) == BCR_f)       //  Condición de extremo izquierdo en la derivada y extremo derecho en la función
    {
        matrixEDP_T diagMA(dim-1,3);
        vectorEDP_T b(dim-1), solV(dim-1);
        
        diagMA(0,0) = 0.0;
        diagMA(0,1) = 2.0 - h*h*B(x(0));
        diagMA(0,2) = -2.0;
        b(0) = -h*h*C(x(0)) + 2.0*h*BCL(x(0),0.0)*(h*A(x(0))/2.0 - 1.0);
        for (int i=1; i<dim-2; i++) {
            tI += h;
            diagMA(i,0) = -1.0 - ((h*A(tI))/2.0);
            diagMA(i,1) = (2.0 - h*h*B(tI));
            diagMA(i,2) = ((h*A(tI))/2.0) - 1.0;
            b(i) = -h*h*C(tI);
        }
        b(dim-2) += (1.0 - (A(x(dim-1)-h)*h)/2.0)*BCR(x(dim-1),0.0);
        
        if ((opt& LUmethod) != 0) {
            solV = diagMA.solveLU3d(b);
        } else if ((opt& GSmethod) != 0) {
            solV = diagMA.solveGS3d(b, err);
        } else {
            solV = diagMA.solveLU3d(b);
        }
        
        sol = y;
        sol.set(0, solV);
    }
    
    return sol;
}



//  -- ECUACIONES DE LAPLACE Y POISSON - SISTEMAS DE 2 DIMENSIONES --
//  Resuelve las ecuaciones diferenciales: ∆U(x,y) = 0 y ∆U(x,y) = F, F ≠ F(U)

//  ECUACIÓN DE LAPLACE
//  Función principal / Método de resolución
matrixEDP_T EDP::solveLAPLACE(unsigned char bc, unsigned char sbc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err, int tol)
{
    int n = cI.rows(), m = cI.columns();
    EDP_T hx = (EDP_T)(x(m-1)-x(0))/(m-1), hy = (EDP_T)(y(n-1)-y(0))/(n-1);
    matrixEDP_T sol(n,m), solOld(n,m);
    
    sol = solOld = cI;
    
    //  Condiciones en los extremos superior e inferior (Condiciones en la función)
    for (int i=0; i<m; i++) {
        if ((bc& BCT_f) != 0) {
            sol(0,i) = BCT(x(i),y(0));
        }
        if ((bc& BCB_f) != 0) {
            sol(n-1,i) = BCB(x(i),y(n-1));
        }
    }
    //  Condiciones en los extremos izquierdo y derecho (Condiciones en la función)
    for (int i=0; i<n; i++) {
        if ((bc& BCL_f) != 0) {
            sol(i,0) = BCL(x(0),y(i));
        }
        if ((bc& BCR_f) != 0) {
            sol(i,m-1) = BCR(x(m-1),y(i));
        }
    }
    
    //  Solución de la ecuación de Poisson
    int ite = 0;
    EDP_T sum, sumOld, relErr;
    vectorEDP_T sup(m-2), supOld(m-2);
    
    //  Solución de la ecuación - Condiciones de contorno normales
    for (int k=0; k<tol; k++) {
        for (int i=1; i<n-1; i++) {
            sum = sumOld = 0.0;
            for (int j=1; j<m-1; j++) {
                //  Condiciones de contorno del borde izquierdo
                if ((bc& BCL_df) != 0)
                    sol(i,0) = (hx*hx*(sol(i+1,0) + sol(i-1,0)) + hy*hy*2*(sol(i,1) + hx*BCL(x(0),y(i))))/(2.0*(hy*hy+hx*hx));
                if ((sbc& BCL_f) != 0)
                    sol(i,0) = SBCL(sol(i,1), sol(i-1,0), sol(i+1,0));
                
                //  Condiciones de contorno del borde superior
                if ((bc& BCT_df) != 0)
                    sol(0,j) = (hx*hx*2*(sol(1,j) + hy*BCT(x(j),y(0))) + hy*hy*(sol(0,j+1) + sol(0,j-1)))/(2.0*(hy*hy+hx*hx));
                if ((sbc& BCT_f) != 0)
                    sol(0,j) = SBCT(sol(0,j-1), sol(0,j+1), sol(1,j));
                
                //  Condiciones de contorno del interior
                if ((sbc& BCI_f) != 0) {
                    sol(i,j) = SBCI(sol(i,j-1), sol(i,j+1), sol(i-1,j), sol(i+1,j));
                } else {
                    sol(i,j) = (hx*hx*(sol(i+1,j) + sol(i-1,j)) + hy*hy*(sol(i,j+1) + sol(i,j-1)))/(2.0*(hy*hy+hx*hx));
                }
                
                //  Condiciones de contorno del borde inferior
                if ((bc& BCB_df) != 0)
                    sol(n-1,j) = (hx*hx*2*(sol(n-2,j) + hy*BCB(x(j),y(n-1))) + hy*hy*(sol(n-1,j+1) + sol(n-1,j-1)))/(2.0*(hy*hy+hx*hx));
                if ((sbc& BCB_f) != 0)
                    sol(n-1,j) = SBCB(sol(n-1,j-1), sol(n-1,j+1), sol(n-2,j));
                
                //  Condiciones de contorno del borde derecho
                if ((bc& BCR_df) != 0)
                    sol(i,m-1) = (hx*hx*(sol(i+1,m-1) + sol(i-1,m-1)) + hy*hy*2*(sol(i,m-2) + hx*BCR(x(m-1),y(i))))/(2.0*(hy*hy+hx*hx));
                if ((sbc& BCR_f) != 0)
                    sol(i,m-1) = SBCR(sol(i,m-2), sol(i-1,m-1), sol(i+1,m-1));
                
                //  Ajuste de las esquinas
                if ((bc& BCL_df) != 0 and (bc& BCT_df) != 0)
                    sol(0,0) = (2.0*sol(0,1)-sol(0,2) + 2.0*sol(1,0)-sol(2,0))/2.0;
                if ((bc& BCT_df) != 0 and (bc& BCR_df) != 0)
                    sol(0,m-1) = (2.0*sol(0,m-2)-sol(0,m-3) + 2.0*sol(1,m-1)-sol(2,m-1))/2.0;
                if ((bc& BCL_df) != 0 and (bc& BCB_df) != 0)
                    sol(n-1,0) = (2.0*sol(n-1,1)-sol(n-1,2) + 2.0*sol(n-2,0)-sol(n-3,0))/2.0;
                if ((bc& BCB_df) != 0 and (bc& BCR_df) != 0)
                    sol(n-1,m-1) = (2.0*sol(n-1,m-2)-sol(n-1,m-3) + 2.0*sol(n-2,m-1)-sol(n-3,m-1))/2.0;
                
                sum += fabs(sol(i,j));
                sumOld += fabs(solOld(i,j));
                solOld(i,j) = sol(i,j);
            }
            
            sup(i-1) = sum;
            supOld(i-1) = sumOld;
            relErr = fabs(sup.max()-supOld.max());
        }
        
        ite++;
        if (relErr < err) {
            break;
        }
    }
    
    if ((opt& ITERATIONS) != 0)
        cout << "\n [Información Laplace]: Se han realizado " << scientific << ite << " iteraciones.\n\n";
    
    if (ite == tol)
        cout << "\n [Información Laplace]: Se ha excedido el número máximo de iteraciones: " << scientific << tol << ", esto quiere decir que el resultado no ha convergido por debajo de la variación entre iteraciones.\n\n";
    
    if ((opt& PLOT_MATLAB) != 0)
        printMATLAB("", "data", DESKTOP, x, y, sol);
    
    return sol;
}

//  Funciones de llamada
matrixEDP_T EDP::solveLAPLACE(unsigned char bc, unsigned char sbc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err)
{
    return solveLAPLACE(bc, sbc, opt, x, y, cI, err, tolDef);
}
matrixEDP_T EDP::solveLAPLACE(unsigned char bc, unsigned char sbc, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err, int tol)
{
    return solveLAPLACE(bc, sbc, 0, x, y, cI, err, tol);
}
matrixEDP_T EDP::solveLAPLACE_NBC(unsigned char bc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err)
{
    return solveLAPLACE(bc, 0, opt, x, y, cI, err, tolDef);
}
matrixEDP_T EDP::solveLAPLACE_SBC(unsigned char sbc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err)
{
    return solveLAPLACE(0, sbc, opt, x, y, cI, err, tolDef);
}
matrixEDP_T EDP::solveLAPLACE_NBC(unsigned char bc, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err)
{
    return solveLAPLACE(bc, 0, 0, x, y, cI, err, tolDef);
}
matrixEDP_T EDP::solveLAPLACE_SBC(unsigned char sbc, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err)
{
    return solveLAPLACE(0, sbc, 0, x, y, cI, err, tolDef);
}


//  ECUACIÓN DE POISSON
//  Función principal / Método de resolución
matrixEDP_T EDP::solvePOISSON(unsigned char bc, unsigned char sbc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err, int tol)
{
    int n = cI.rows(), m = cI.columns();
    EDP_T hx = (EDP_T)(x(m-1)-x(0))/(m-1), hy = (EDP_T)(y(n-1)-y(0))/(n-1);
    matrixEDP_T sol(n,m), solOld(n,m);
    
    sol = solOld = cI;
    
    //  Condiciones en los extremos superior e inferior (Condiciones en la función)
    for (int i=0; i<m; i++) {
        if ((bc& BCT_f) != 0) {
            sol(0,i) = BCT(x(i),y(0));
        }
        if ((bc& BCB_f) != 0) {
            sol(n-1,i) = BCB(x(i),y(n-1));
        }
    }
    //  Condiciones en los extremos izquierdo y derecho (Condiciones en la función)
    for (int i=0; i<n; i++) {
        if ((bc& BCL_f) != 0) {
            sol(i,0) = BCL(x(0),y(i));
        }
        if ((bc& BCR_f) != 0) {
            sol(i,m-1) = BCR(x(m-1),y(i));
        }
    }
    
    //  Solución de la ecuación de Poisson
    int ite = 0;
    EDP_T sum, sumOld, relErr;
    vectorEDP_T sup(m-2), supOld(m-2);
    
    //  Solución de la ecuación - Condiciones de contorno normales
    for (int k=0; k<tol; k++) {
        for (int i=1; i<n-1; i++) {
            sum = sumOld = 0.0;
            for (int j=1; j<m-1; j++) {
                //  Condiciones de contorno del borde izquierdo
                if ((bc& BCL_df) != 0)
                    sol(i,0) = (hx*hx*(sol(i+1,0) + sol(i-1,0)) + hy*hy*2*(sol(i,1) + hx*BCL(x(0),y(i))) - hy*hy*hx*hx*F(x(0),y(i)))/(2.0*(hy*hy+hx*hx));
                if ((sbc& BCL_f) != 0)
                    sol(i,0) = SBCL(sol(i,1), sol(i-1,0), sol(i+1,0));
                
                //  Condiciones de contorno del borde superior
                if ((bc& BCT_df) != 0)
                    sol(0,j) = (hx*hx*2*(sol(1,j) + hy*BCT(x(j),y(0))) + hy*hy*(sol(0,j+1) + sol(0,j-1)) - hy*hy*hx*hx*F(x(j),y(0)))/(2.0*(hy*hy+hx*hx));
                if ((sbc& BCT_f) != 0)
                    sol(0,j) = SBCT(sol(0,j-1), sol(0,j+1), sol(1,j));
                
                //  Condiciones de contorno del interior
                if ((sbc& BCI_f) != 0) {
                    sol(i,j) = SBCI(sol(i,j-1), sol(i,j+1), sol(i-1,j), sol(i+1,j));
                } else {
                    sol(i,j) = (hx*hx*(sol(i+1,j) + sol(i-1,j)) + hy*hy*(sol(i,j+1) + sol(i,j-1)) - hy*hy*hx*hx*F(x(j),y(i)))/(2.0*(hy*hy+hx*hx));
                }
                
                //  Condiciones de contorno del borde inferior
                if ((bc& BCB_df) != 0)
                    sol(n-1,j) = (hx*hx*2*(sol(n-2,j) + hy*BCB(x(j),y(n-1))) + hy*hy*(sol(n-1,j+1) + sol(n-1,j-1)) - hy*hy*hx*hx*F(x(j),y(n-1)))/(2.0*(hy*hy+hx*hx));
                if ((sbc& BCB_f) != 0)
                    sol(n-1,j) = SBCB(sol(n-1,j-1), sol(n-1,j+1), sol(n-2,j));
                
                //  Condiciones de contorno del borde derecho
                if ((bc& BCR_df) != 0)
                    sol(i,m-1) = (hx*hx*(sol(i+1,m-1) + sol(i-1,m-1)) + hy*hy*2*(sol(i,m-2) + hx*BCR(x(m-1),y(i))) - hy*hy*hx*hx*F(x(m-1),y(i)))/(2.0*(hy*hy+hx*hx));
                if ((sbc& BCR_f) != 0)
                    sol(i,m-1) = SBCR(sol(i,m-2), sol(i-1,m-1), sol(i+1,m-1));
                
                //  Ajuste de las esquinas
                if ((bc& BCL_df) != 0 and (bc& BCT_df) != 0)
                    sol(0,0) = (2.0*sol(0,1)-sol(0,2) + 2.0*sol(1,0)-sol(2,0))/2.0;
                if ((bc& BCT_df) != 0 and (bc& BCR_df) != 0)
                    sol(0,m-1) = (2.0*sol(0,m-2)-sol(0,m-3) + 2.0*sol(1,m-1)-sol(2,m-1))/2.0;
                if ((bc& BCL_df) != 0 and (bc& BCB_df) != 0)
                    sol(n-1,0) = (2.0*sol(n-1,1)-sol(n-1,2) + 2.0*sol(n-2,0)-sol(n-3,0))/2.0;
                if ((bc& BCB_df) != 0 and (bc& BCR_df) != 0)
                    sol(n-1,m-1) = (2.0*sol(n-1,m-2)-sol(n-1,m-3) + 2.0*sol(n-2,m-1)-sol(n-3,m-1))/2.0;
                
                sum += fabs(sol(i,j));
                sumOld += fabs(solOld(i,j));
                solOld(i,j) = sol(i,j);
            }
            
            sup(i-1) = sum;
            supOld(i-1) = sumOld;
            relErr = fabs(sup.max()-supOld.max());
        }
        
        ite++;
        if (relErr < err) {
            break;
        }
    }
    
    if ((opt& ITERATIONS) != 0)
        cout << "\n [Información Poisson]: Se han realizado " << scientific << ite << " iteraciones.\n\n";
    
    if (ite == tol)
        cout << "\n [Información Poisson]: Se ha excedido el número máximo de iteraciones: " << scientific << tol << ", esto quiere decir que el resultado no ha convergido por debajo de la variación entre iteraciones.\n\n";
    
    if ((opt& PLOT_MATLAB) != 0)
        printMATLAB("", "data", DESKTOP, x, y, sol);
    
    if ((opt& SAVE_DATA) != 0)
        saveDATA("", "data", DESKTOP, x, y, sol);
    
    return sol;
}

//  Funciones de llamada
matrixEDP_T EDP::solvePOISSON(unsigned char bc, unsigned char sbc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err)
{
    return solvePOISSON(bc, sbc, opt, x, y, cI, err, tolDef);
}
matrixEDP_T EDP::solvePOISSON(unsigned char bc, unsigned char sbc, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err, int tol)
{
    return solvePOISSON(bc, sbc, 0, x, y, cI, err, tol);
}
matrixEDP_T EDP::solvePOISSON_NBC(unsigned char bc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err)
{
    return solvePOISSON(bc, 0, opt, x, y, cI, err, tolDef);
}
matrixEDP_T EDP::solvePOISSON_SBC(unsigned char sbc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err)
{
    return solvePOISSON(0, sbc, opt, x, y, cI, err, tolDef);
}
matrixEDP_T EDP::solvePOISSON_NBC(unsigned char bc, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err)
{
    return solvePOISSON(bc, 0, 0, x, y, cI, err, tolDef);
}
matrixEDP_T EDP::solvePOISSON_SBC(unsigned char sbc, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI, EDP_T err)
{
    return solvePOISSON(0, sbc, 0, x, y, cI, err, tolDef);
}



//  -- ECUACIONES DE ONDA Y DE CALOR --
//  Resuelve las EDPs de onda: ∂²u/∂t² = Q()·∆u y de calor: ∂u/∂t = Q()·∆u

//  ECUACIÓN DE ONDA
//  1 Dimensión
vectorEDP_T EDP::solveWAVE(unsigned char bc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& cI, vectorEDP_T& cId)
{
    int n = (int)x.dim();
    vectorEDP_T sol(n);
    EDP_T dx = (EDP_T)abs((EDP_T)(x(n-1) - x(0))/(n-1));
    
    if (dt > dx/sqrt(Q1D(x(0)))) {
        dt = dx/sqrt(Q1D(x(0))) * 0.9;
        cout << EDPwarning << "solveWAVE(bc, opt, x, cI, cId)] - El diferencial de tiempo era demasiado grande para obtener buenos resultados, se ha cambiado por: " << dt << endl;
    }
    
    EDP_T dtx = ((dt*dt)/(dx*dx));
    sol = cI;
    
    if (!initEDP) {
        initEDP = true;
        
        old1D = old1D.zero(n);
        if ((bc& BCL_df) != 0)
            sol(0) = cId(0)*dt + (1.0 - Q1D(x(0))*dtx)*cI(0) + Q1D(x(0))*dtx*(cI(1) - dx*BCL(x(0), 0.0));
        
        for (int i=1; i<n-1; i++) {
            sol(i) = cId(i)*dt + (1.0 - Q1D(x(i))*dtx)*cI(i) + (Q1D(x(i))*dtx/2.0)*(cI(i+1) + cI(i-1));
        }
        
        if ((bc& BCR_df) != 0)
            sol(n-1) = cId(n-1)*dt + (1.0 - Q1D(x(n-1))*dtx)*cI(n-1) + Q1D(x(n-1))*dtx*(cI(n-2) + dx*BCR(x(n-1), 0.0));
    } else {
        if ((bc& BCL_df) != 0)
            sol(0) = 2.0*(1.0 - Q1D(x(0))*dtx)*cI(0) + 2.0*Q1D(x(0))*dtx*(cI(1) - dx*BCL(x(0), 0.0)) - old1D(0);
        
        for (int i=1; i<n-1; i++)
            sol(i) = 2.0*(1.0 - Q1D(x(i))*dtx)*cI(i) + Q1D(x(i))*dtx*(cI(i+1) + cI(i-1)) - old1D(i);
        
        if ((bc& BCR_df) != 0)
            sol(n-1) = 2.0*(1.0 - Q1D(x(n-1))*dtx)*cI(n-1) + 2.0*Q1D(x(n-1))*dtx*(cI(n-2) + dx*BCR(x(n-1), 0.0)) - old1D(n-1);
    }
    
    old1D = cI;
    time += dt;
    
    return sol;
}

vectorEDP_T EDP::solveWAVE(unsigned char bc, vectorEDP_T& x, vectorEDP_T& cI, vectorEDP_T& cId)
{
    return solveWAVE(bc, 0, x, cI, cId);
}

//  2 Dimensiones
matrixEDP_T EDP::solveWAVE(unsigned char bc, unsigned char opt, vectorEDP_T &x, vectorEDP_T &y, matrixEDP_T &cI, matrixEDP_T &cId, matrix<bool> &fixed)
{
    int n = (int)y.dim();
    int m = (int)x.dim();
    matrixEDP_T sol(n,m);
    EDP_T dx = (EDP_T)abs((EDP_T)(x(m-1) - x(0))/(m-1));
    EDP_T dy = (EDP_T)abs((EDP_T)(y(n-1) - y(0))/(n-1));
    
    if (dt > dx*dy/(sqrt(Q2D(x(0),y(0))*(dx*dx+dy*dy)))) {
        dt = dx*dy/(sqrt(Q2D(x(0),y(0))*(dx*dx+dy*dy))) * 0.9;
        cout << EDPwarning << "solveWAVE(bc, opt, x, cI, cId)] - El diferencial de tiempo era demasiado grande para obtener buenos resultados, se ha cambiado por: " << dt << endl;
    }
    
    sol = cI;
    
    if (!initEDP) {
        initEDP = true;
        
        old2D = old2D.zero(n,m);
        if ((bc& BCT_df) != 0)  //  Condición en el borde superior de la membrana
            for (int j=1; j<m-1; j++) {
                if (!fixed(0,j))
                    sol(0,j) = cI(0,j) + dt*cId(0,j) + Q2D(x(j),y(0))*dt*dt/2.0*((2.0*cI(1,j)-2.0*dy*BCT(x(j),y(0))-2.0*cI(0,j))/(dy*dy) + (cI(0,j+1)-2.0*cI(0,j)+cI(0,j-1))/(dx*dx));
            }
        
        if ((bc& BCL_df) != 0)//  Condición en el borde izquierdo de la membrana
            for (int i=1; i<n-1; i++) {
                if (!fixed(i,0))
                    sol(i,0) = cI(i,0) + dt*cId(i,0) + Q2D(x(0),y(i))*dt*dt/2.0*((cI(i+1,0)-2.0*cI(i,0)+cI(i-1,0))/(dy*dy) + (2.0*cI(i,1)-2.0*dx*BCL(x(0),y(i))-2.0*cI(i,0))/(dx*dx));
            }
        
        for (int i=1; i<n-1; i++) {
            for (int j=1; j<m-1; j++) {
                if (!fixed(i,j))
                    sol(i,j) = cI(i,j) + dt*cId(i,j) + Q2D(x(j),y(i))*dt*dt/2.0*((cI(i+1,j)-2.0*cI(i,j)+cI(i-1,j))/(dy*dy) + (cI(i,j+1)-2.0*cI(i,j)+cI(i,j-1))/(dx*dx));
            }
        }
        
        if ((bc& BCR_df) != 0)  //  Condición en el borde derecho de la membrana
            for (int i=1; i<n-1; i++) {
                if (!fixed(i,m-1))
                    sol(i,m-1) = cI(i,m-1) + dt*cId(i,m-1) + Q2D(x(m-1),y(i))*dt*dt/2.0*((cI(i+1,m-1)-2.0*cI(i,m-1)+cI(i-1,m-1))/(dy*dy) + (2.0*cI(i,m-2)+2.0*dx*BCR(x(m-1),y(i))-2.0*cI(i,m-1))/(dx*dx));
            }
        
        if ((bc& BCB_df) != 0)  //  Condición en el borde inferior de la membrana
            for (int j=1; j<m-1; j++) {
                if (!fixed(n-1,j))
                    sol(n-1,j) = cI(n-1,j) + dt*cId(n-1,j) + Q2D(x(j),y(n-1))*dt*dt/2.0*((2.0*cI(n-2,j)+2.0*dy*BCB(x(j),y(n-1))-2.0*cI(n-1,j))/(dy*dy) + (cI(n-1,j+1)-2.0*cI(n-1,j)+cI(n-1,j-1))/(dx*dx));
            }
        
        if ((bc& BCL_df) != 0 and (bc& BCT_df) != 0 and !fixed(0,0))
            sol(0,0) = cI(0,0) + dt*cId(0,0) + Q2D(x(0),y(0))*dt*dt*((cI(1,0)-cI(0,0)-dy*BCT(x(0),y(0)))/(dy*dy) + (cI(0,1)-cI(0,0)-dx*BCL(x(0),y(0)))/(dx*dx));
                
        if ((bc& BCT_df) != 0 and (bc& BCR_df) != 0 and !fixed(0,m-1))
            sol(0,m-1) = cI(0,m-1) + dt*cId(0,m-1) + Q2D(x(m-1),y(0))*dt*dt*((cI(1,m-1)-cI(0,m-1)-dy*BCT(x(m-1),y(0)))/(dy*dy) + (cI(0,m-2)-cI(0,m-1)+dx*BCR(x(m-1),y(0)))/(dx*dx));
        
        if ((bc& BCL_df) != 0 and (bc& BCB_df) != 0 and !fixed(n-1,0))
            sol(n-1,0) = cI(n-1,0) + dt*cId(n-1,0) + Q2D(x(0),y(n-1))*dt*dt*((cI(n-2,0)-cI(n-1,0)+dy*BCB(x(0),y(n-1)))/(dy*dy) + (cI(n-1,1)-cI(n-1,0)-dx*BCL(x(0),y(n-1)))/(dx*dx));
        
        if ((bc& BCB_df) != 0 and (bc& BCR_df) != 0 and !fixed(n-1,m-1))
            sol(n-1,m-1) = cI(n-1,m-1) + dt*cId(n-1,m-1) + Q2D(x(m-1),y(n-1))*dt*dt*((cI(n-2,m-1)-cI(n-1,m-1)+dy*BCB(x(m-1),y(n-1)))/(dy*dy) + (cI(n-1,m-2)-cI(n-1,m-1)+dx*BCR(x(m-1),y(n-1)))/(dx*dx));
        
    } else {
        if ((bc& BCT_df) != 0)  //  Condición en el borde superior de la membrana
            for (int j=1; j<m-1; j++) {
                if (!fixed(0,j))
                    sol(0,j) = 2.0*cI(0,j) + dt*dt*Q2D(x(j),y(0))*((2.0*cI(1,j)-2.0*dy*BCT(x(j),y(0))-2.0*cI(0,j))/(dy*dy) + (cI(0,j+1)-2.0*cI(0,j)+cI(0,j-1))/(dx*dx)) - old2D(0,j);
            }
        
        if ((bc& BCL_df) != 0)  //  Condición en el borde izquierdo de la membrana
            for (int i=1; i<n-1; i++) {
                if (!fixed(i,0))
                    sol(i,0) = 2.0*cI(i,0) + dt*dt*Q2D(x(0),y(i))*((cI(i+1,0)-2.0*cI(i,0)+cI(i-1,0))/(dy*dy) + (2.0*cI(i,1)-2.0*dx*BCL(x(0),y(i))-2.0*cI(i,0))/(dx*dx)) - old2D(i,0);
            }
        
        
        for (int i=1; i<n-1; i++) {
            for (int j=1; j<m-1; j++) {
                if (!fixed(i,j))
                    sol(i,j) = 2.0*cI(i,j) + dt*dt*Q2D(x(j),y(i))*((cI(i+1,j)-2.0*cI(i,j)+cI(i-1,j))/(dy*dy) + (cI(i,j+1)-2.0*cI(i,j)+cI(i,j-1))/(dx*dx)) - old2D(i,j);
            }
        }
        
        
        if ((bc& BCR_df) != 0)  //  Condición en el borde derecho de la membrana
            for (int i=1; i<n-1; i++) {
                if (!fixed(i,m-1))
                    sol(i,m-1) = 2.0*cI(i,m-1) + dt*dt*Q2D(x(m-1),y(i))*((cI(i+1,m-1)-2.0*cI(i,m-1)+cI(i-1,m-1))/(dy*dy) + (2.0*cI(i,m-2)+2.0*dx*BCR(x(m-1),y(i))-2.0*cI(i,m-1))/(dx*dx)) - old2D(i,m-1);
            }
        
        if ((bc& BCB_df) != 0)  //  Condición en el borde inferior de la membrana
            for (int j=1; j<m-1; j++) {
                if (!fixed(n-1,j))
                    sol(n-1,j) = 2.0*cI(n-1,j) + dt*dt*Q2D(x(j),y(n-1))*((2.0*cI(n-2,j)+2.0*dy*BCB(x(j),y(n-1))-2.0*cI(n-1,j))/(dy*dy) + (cI(n-1,j+1)-2.0*cI(n-1,j)+cI(n-1,j-1))/(dx*dx)) - old2D(n-1,j);
            }
        
        if ((bc& BCL_df) != 0 and (bc& BCT_df) != 0 and !fixed(0,0))
            sol(0,0) = 2.0*(cI(0,0) + dt*dt*Q2D(x(0),y(0))*((cI(1,0)-cI(0,0)-dy*BCT(x(0),y(0)))/(dy*dy) + (cI(0,1)-cI(0,0)-dx*BCL(x(0),y(0)))/(dx*dx))) - old2D(0,0);
        
        if ((bc& BCT_df) != 0 and (bc& BCR_df) != 0 and !fixed(0,m-1))
            sol(0,m-1) = 2.0*(cI(0,m-1) + dt*dt*Q2D(x(m-1),y(0))*((cI(1,m-1)-cI(0,m-1)-dy*BCT(x(m-1),y(0)))/(dy*dy) + (cI(0,m-2)-cI(0,m-1)+dx*BCR(x(m-1),y(0)))/(dx*dx))) - old2D(0,m-1);
        
        if ((bc& BCL_df) != 0 and (bc& BCB_df) != 0 and !fixed(n-1,0))
            sol(n-1,0) = 2.0*(cI(n-1,0) + dt*dt*Q2D(x(0),y(n-1))*((cI(n-2,0)-cI(n-1,0)+dy*BCB(x(0),y(n-1)))/(dy*dy) + (cI(n-1,1)-cI(n-1,0)-dx*BCL(x(0),y(n-1)))/(dx*dx))) - old2D(n-1,0);
        
        if ((bc& BCB_df) != 0 and (bc& BCR_df) != 0 and !fixed(n-1,m-1))
            sol(n-1,m-1) = 2.0*(cI(n-1,m-1) + dt*dt*Q2D(x(m-1),y(n-1))*((cI(n-2,m-1)-cI(n-1,m-1)+dy*BCB(x(m-1),y(n-1)))/(dy*dy) + (cI(n-1,m-2)-cI(n-1,m-1)+dx*BCR(x(m-1),y(n-1)))/(dx*dx))) - old2D(n-1,m-1);
    }
    
    old2D = cI;
    time += (EDP_T)dt;
    
    if ((opt& SAVE_DATA) != 0) {
        if (pathEDP != "")
            saveDATA(pathEDP, "WAVE_Equation.csv", x, y, sol);
        else 
            saveDATA("WAVE_Equation.csv", DESKTOP, x, y, sol);
    }
    
    return sol;
}

matrixEDP_T EDP::solveWave(unsigned char bc, vectorEDP_T &x, vectorEDP_T &y, matrixEDP_T &cI, matrixEDP_T &cId, matrix<bool> &fixed)
{
    return solveWAVE(bc, 0, x, y, cI, cId, fixed);
}

matrixEDP_T EDP::solveWAVE(unsigned char bc, vectorEDP_T &x, vectorEDP_T &y, matrixEDP_T &cI, matrixEDP_T &cId)
{
    if (!initEDP)
        fixedEDP = zero<bool>(y.dim(), x.dim());
    
    return solveWAVE(bc, 0, x, y, cI, cId, fixedEDP);
}


//  ECUACIÓN DEL CALOR
//  1 Dimensión
vectorEDP_T EDP::solveHEAT(unsigned char bc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y)
{
    int n = (int)x.dim();
    vectorEDP_T sol(n);
    
    if ((bc& BCL_f) != 0 and (bc& BCR_f) != 0) {
        matrixEDP_T diagLU(n-2,3);
        vectorEDP_T b(n-2), solV(n-2);
        EDP_T dx = (EDP_T)abs((EDP_T)(x(n-1) - y(0))/(n-1));
        EDP_T dtx = dt/(dx*dx);
        
        if (theta < 0 || theta > 1) {
            theta = (EDP_T)2/3;
        }
        
        for (int i=0; i<n-2; i++) {
            diagLU(i,0) = -Q1D(x(i))*dtx*theta;
            diagLU(i,1) = 1.0 + 2.0*Q1D(x(i))*dtx*theta;
            diagLU(i,2) = -Q1D(x(i))*dtx*theta;
            b(i) = Q1D(x(i))*dtx*(1.0-theta)*y(i+2) + (1.0 - 2.0*Q1D(x(i))*dtx*(1.0-theta))*y(i+1) + Q1D(x(i))*dtx*(1.0-theta)*y(i);
        }
        
        solV = diagLU.solveLU3d(b);
        
        sol = y;
        sol.set(1, solV);
    }
    
    if ((bc& BCL_f) != 0 and (bc& BCR_df) != 0) {
        matrixEDP_T diagLU(n-1,3);
        vectorEDP_T b(n-1), solV(n-1);
        EDP_T dx = (EDP_T)abs((EDP_T)(y(n-1) - y(0))/(n-1));
        EDP_T dtx = dt/(dx*dx);
        
        if (theta < 0 || theta > 1) {
            theta = (EDP_T)2/3;
        }
        
        for (int i=0; i<n-2; i++) {
            diagLU(i,0) = -Q1D(x(i))*dtx*theta;
            diagLU(i,1) = 1.0 + 2.0*Q1D(x(i))*dtx*theta;
            diagLU(i,2) = -Q1D(x(i))*dtx*theta;
            b(i) = Q1D(x(i))*dtx*(1.0-theta)*y(i+2) + (1.0 - 2.0*Q1D(x(i))*dtx*(1.0-theta))*y(i+1) + Q1D(x(i))*dtx*(1.0-theta)*y(i);
        }
        
        diagLU(n-2,0) = -2*Q1D(x(n-1))*dtx*theta;
        diagLU(n-2,1) = 1.0 + 2.0*Q1D(x(n-1))*dtx*theta;
        diagLU(n-2,2) = 0.0;
        b(n-2) = 2.0*Q1D(x(n-1))*dtx*(1.0-theta)*y(n-2) + (1.0 - 2.0*Q1D(x(n-1))*dtx*(1.0-theta))*y(n-1) + 2.0*Q1D(x(n-1))*dtx*(1.0-2.0*theta)*dx*BCR(x(n-1),0.0);
        
        solV = diagLU.solveLU3d(b);
        
        sol = y;
        y.set(1, solV);
    }
    
    if ((bc& BCL_df) != 0 and (bc& BCR_f) != 0) {
        matrixEDP_T diagLU(n-1,3);
        vectorEDP_T b(n-1), solV(n-1);
        EDP_T dx = (EDP_T)abs((EDP_T)(x(n-1) - x(0))/(n-1));
        EDP_T dtx = dt/(dx*dx);
        
        if (theta < 0 || theta > 1) {
            theta = (EDP_T)2/3;
        }
        
        diagLU(0,0) = 0.0;
        diagLU(0,1) = 1.0 + 2.0*Q1D(x(0))*dtx*theta;
        diagLU(0,2) = -2.0*Q1D(x(0))*dtx*theta;
        b(0) = 2.0*Q1D(x(0))*dtx*(1.0-theta)*y(1) + (1.0 - 2.0*Q1D(x(0))*dtx*(1.0-theta))*y(0) - 2.0*Q1D(x(0))*dtx*(1.0-2.0*theta)*dx*BCL(x(0),0.0);
        
        for (int i=1; i<n-1; i++) {
            diagLU(i,0) = -Q1D(x(i))*dtx*theta;
            diagLU(i,1) = 1.0 + 2.0*Q1D(x(i))*dtx*theta;
            diagLU(i,2) = -Q1D(x(i))*dtx*theta;
            b(i) = Q1D(x(i))*dtx*(1.0-theta)*y(i+1) + (1.0 - 2.0*Q1D(x(i))*dtx*(1.0-theta))*y(i) + Q1D(x(i))*dtx*(1.0-theta)*y(i-1);
        }
        
        solV = diagLU.solveLU3d(b);
        
        sol = y;
        sol.set(0, solV);
    }
    
    if ((bc& BCL_df) != 0 and (bc& BCR_df) != 0) {
        matrixEDP_T diagLU(n,3);
        vectorEDP_T b(n), solV(n);
        EDP_T dx = (EDP_T)abs((EDP_T)(x(n-1) - x(0))/(n-1));
        EDP_T dtx = dt/(dx*dx);
        
        if (theta < 0 || theta > 1) {
            theta = (EDP_T)2/3;
        }
        
        diagLU(0,0) = 0.0;
        diagLU(0,1) = 1.0 + 2.0*Q1D(x(0))*dtx*theta;
        diagLU(0,2) = -2.0*Q1D(x(0))*dtx*theta;
        b(0) = 2.0*Q1D(x(0))*dtx*(1.0-theta)*y(1) + (1.0 - 2.0*Q1D(x(0))*dtx*(1.0-theta))*y(0) - 2.0*Q1D(x(0))*dtx*(1.0-2.0*theta)*dx*BCL(x(0),0.0);
        
        for (int i=1; i<n-1; i++) {
            diagLU(i,0) = -Q1D(x(i))*dtx*theta;
            diagLU(i,1) = 1.0 + 2.0*Q1D(x(i))*dtx*theta;
            diagLU(i,2) = -Q1D(x(i))*dtx*theta;
            b(i) = Q1D(x(i))*dtx*(1.0-theta)*y(i+1) + (1.0 - 2.0*Q1D(x(i))*dtx*(1.0-theta))*y(i) + Q1D(x(i))*dtx*(1.0-theta)*y(i-1);
        }
        
        diagLU(n-1,0) = -2*Q1D(x(n-1))*dtx*theta;
        diagLU(n-1,1) = 1.0 + 2.0*Q1D(x(n-1))*dtx*theta;
        diagLU(n-1,2) = 0.0;
        b(n-1) = 2.0*Q1D(x(n-1))*dtx*(1.0-theta)*y(n-2) + (1.0 - 2.0*Q1D(x(n-1))*dtx*(1.0-theta))*y(n-1) + 2.0*Q1D(x(n-1))*dtx*(1.0-2.0*theta)*dx*BCR(x(n-1),0.0);
        
        solV = diagLU.solveLU3d(b);
        
        sol = solV;
    }
    
    time += dt;
    
    return sol;
}

vectorEDP_T EDP::solveHEAT(unsigned char bc, vectorEDP_T& x, vectorEDP_T& y)
{
    return solveHEAT(bc, 0, x, y);
}

//  2 Dimensiones
matrixEDP_T EDP::solveHEAT(unsigned char bc, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI)
{
    int n = (int)y.dim();
    int m = (int)x.dim();
    matrixEDP_T sol(n,m);
    EDP_T dx = (EDP_T)abs((EDP_T)(x(m-1) - x(0))/(m-1));
    EDP_T dy = (EDP_T)abs((EDP_T)(y(n-1) - y(0))/(n-1));
    
    if (dt > dx*dx*dy*dy/(dx*dx+dy*dy)*1/(4*Q2D(x(0),y(0)))) {
        dt = dx*dx*dy*dy/(dx*dx+dy*dy)*1/(4*Q2D(x(0),y(0))) * 0.9;
        cout << " El diferencial de tiempo era demasiado grande para obtener buenos resultados, se ha cambiado por: " << dt << endl;
    }
    
    sol = cI;
    
    if ((bc& BCT_df) != 0)  //  Condición en el borde superior de la membrana
        for (int j=1; j<m-1; j++) {
            sol(0,j) = cI(0,j) + dt*Q2D(x(j),y(0))*((2.0*cI(1,j)-2.0*dy*BCT(x(j),y(0))-2.0*cI(0,j))/(dy*dy) + (cI(0,j+1)-2.0*cI(0,j)+cI(0,j-1))/(dx*dx));
        }
    
    if ((bc& BCL_df) != 0)  //  Condición en el borde izquierdo de la membrana
        for (int i=1; i<n-1; i++) {
            sol(i,0) = cI(i,0) + dt*Q2D(x(0),y(i))*((cI(i+1,0)-2.0*cI(i,0)+cI(i-1,0))/(dy*dy) + (2.0*cI(i,1)-2.0*dx*BCL(x(0),y(i))-2.0*cI(i,0))/(dx*dx));
        }
    
    for (int i=1; i<n-1; i++) {
        for (int j=1; j<m-1; j++) {
            sol(i,j) = cI(i,j) + dt*Q2D(x(j),y(i))*((cI(i+1,j)-2.0*cI(i,j)+cI(i-1,j))/(dy*dy) + (cI(i,j+1)-2.0*cI(i,j)+cI(i,j-1))/(dx*dx));
        }
    }
    
    if ((bc& BCR_df) != 0)  //  Condición en el borde derecho de la membrana
        for (int i=1; i<n-1; i++) {
            sol(i,m-1) = cI(i,m-1) + dt*Q2D(x(m-1),y(i))*((cI(i+1,m-1)-2.0*cI(i,m-1)+cI(i-1,m-1))/(dy*dy) + (2.0*cI(i,m-2)+2.0*dx*BCR(x(m-1),y(i))-2.0*cI(i,m-1))/(dx*dx));
        }
    
    if ((bc& BCB_df) != 0)  //  Condición en el borde inferior de la membrana
        for (int j=1; j<m-1; j++) {
            sol(n-1,j) = cI(n-1,j) + dt*Q2D(x(j),y(n-1))*((2.0*cI(n-2,j)+2.0*dy*BCB(x(j),y(n-1))-2.0*cI(n-1,j))/(dy*dy) + (cI(n-1,j+1)-2.0*cI(n-1,j)+cI(n-1,j-1))/(dx*dx));
        }
    
    if ((bc& BCL_df) != 0 and (bc& BCT_df) != 0)
        sol(0,0) = (2.0*sol(0,1)-sol(0,2) + 2.0*sol(1,0)-sol(2,0))/2.0;
    if ((bc& BCT_df) != 0 and (bc& BCR_df) != 0)
        sol(0,m-1) = (2.0*sol(0,m-2)-sol(0,m-3) + 2.0*sol(1,m-1)-sol(2,m-1))/2.0;
    if ((bc& BCL_df) != 0 and (bc& BCB_df) != 0)
        sol(n-1,0) = (2.0*sol(n-1,1)-sol(n-1,2) + 2.0*sol(n-2,0)-sol(n-3,0))/2.0;
    if ((bc& BCB_df) != 0 and (bc& BCR_df) != 0)
        sol(n-1,m-1) = (2.0*sol(n-1,m-2)-sol(n-1,m-3) + 2.0*sol(n-2,m-1)-sol(n-3,m-1))/2.0;
    
    time += dt;
    
    if ((opt& SAVE_DATA) != 0) {
        if (pathEDP != "")
            saveDATA(pathEDP, "HEATEquation.csv", x, y, sol);
        else 
            saveDATA("HEATEquation.csv", DESKTOP, x, y, sol);
    }
    
    return sol;
}

matrixEDP_T EDP::solveHEAT(unsigned char bc, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& cI)
{
    return solveHEAT(bc, 0, x, y, cI);
}



//  -- ECUACIONES DEL TIPO -> ∂²u/∂t² = -k²u --
//  Calcula los autovalores y autovectores de dicha ecuación.
vectorEDP_T EDP::eigenVAL_VEC(vectorEDP_T& x, int mode, unsigned char bc, unsigned char opt)
{
    int dim = x.dim();
    EDP_T h = (x(dim-1) - x(0))/(dim-1), length= x(dim-1)-x(0);
    vectorEDP_T sol = zero<EDP_T>(dim), solV, eigVal;
    matrixEDP_T A;
    
    //  Para guardar y cargar datos
    string inPath, fileName;
    if (pathEDP != "")
        inPath = pathEDP;
    else {
        const string home = getenv("HOME");
        inPath = home + "/Desktop/";
    }
    
    //  Posibles casos que se pueden plantear según las condiciones de contorno
    if ((((bc &BCL_df) and (bc &BCR_f)) != 0) || (((bc &BCB_df) and (bc &BCT_f)) != 0)) {
        solV = zero<EDP_T>(dim-1);
        A = zero<EDP_T>(dim-1,dim-1);
        
        A(0,0) = 2.0;
        A(0,1) = -2.0;
        for (int i=1; i<dim-2; i++) {
            A(i,i) = 2.0;
            A(i,i-1) = A(i,i+1) = -1.0;
        }
        A(dim-2,dim-2) = 2.0;
        A(dim-2,dim-3) = -1.0;
        
        eigVal = zero<EDP_T>(dim-1);
        
        if (((bc &BCL_df) and (bc &BCR_f)) != 0)
            fileName = "EV_BCLdf_BCRf.csv";
        else
            fileName = "EV_BCBdf_BCTf.csv";
        
    } else if ((((bc &BCL_f) and (bc &BCR_df)) != 0) || (((bc &BCB_f) and (bc &BCT_df)) != 0)) {
        solV = zero<EDP_T>(dim-1);
        A = zero<EDP_T>(dim-1,dim-1);
        
        A(0,0) = 2.0;
        A(0,1) = -1.0;
        for (int i=1; i<dim-2; i++) {
            A(i,i) = 2.0;
            A(i,i-1) = A(i,i+1) = -1.0;
        }
        A(dim-2,dim-2) = 2.0;
        A(dim-2,dim-3) = -2.0;
        
        eigVal = zero<EDP_T>(dim-1);
        
        if (((bc &BCL_f) and (bc &BCR_df)) != 0)
            fileName = "EV_BCLf_BCRdf.csv";
        else
            fileName = "EV_BCBf_BCTdf.csv";
        
    } else if ((((bc &BCL_df) and (bc &BCR_df)) != 0) || (((bc &BCB_df) and (bc &BCT_df)) != 0)) {
        solV = zero<EDP_T>(dim);
        A = zero<EDP_T>(dim,dim);
        
        A(0,0) = 2.0;
        A(0,1) = -2.0;
        for (int i=1; i<dim-1; i++) {
            A(i,i) = 2.0;
            A(i,i-1) = A(i,i+1) = -1.0;
        }
        A(dim-1,dim-1) = 2.0;
        A(dim-1,dim-2) = -2.0;
        
        eigVal = zero<EDP_T>(dim);
        
        if (((bc &BCL_df) and (bc &BCR_df)) != 0)
            fileName = "EV_BCLdf_BCRdf.csv";
        else
            fileName = "EV_BCBdf_BCTdf.csv";
        
    } else {
        solV = zero<EDP_T>(dim-2);
        A = zero<EDP_T>(dim-2,dim-2);
        
        A(0,0) = 2.0;
        A(0,1) = -1.0;
        for (int i=1; i<dim-3; i++) {
            A(i,i) = 2.0;
            A(i,i-1) = A(i,i+1) = -1.0;
        }
        A(dim-3,dim-3) = 2.0;
        A(dim-3,dim-4) = -1.0;
        
        eigVal = zero<EDP_T>(dim-2);
        
        if (((bc &BCL_f) and (bc &BCR_f)) != 0)
            fileName = "EV_BCLf_BCRf.csv";
        else
            fileName = "EV_BCBf_BCTf.csv";
    }
    
    if ((opt& SAVE_DATA) != 0) {
        cout << "\tCalculando y guardando autovalores... ";
        eigVal = A.eigenValues(20, opt);
        eigVal = eigVal.organize();
        if ((((bc &BCL_df) and (bc &BCR_df)) != 0) || (((bc &BCB_df) and (bc &BCT_df)) != 0))
            eigVal(0) = 0.0;
        const string path = inPath + fileName;
        ofstream out(path.data());
        out.precision(15);
        out << sqrt(eigVal)*length/h;
        cout << "Terminado.\n";
        cout << "\tLos datos se han guardado en: " << path << endl;
    } else if ((opt& IMPORT_DATA) != 0) {
        cout << "\tImportando autovalores... ";
        const string path = inPath + fileName;
        ifstream in(path.data());
        if (in.fail()) {
            cout << "\n\tEl fichero no existe, se van a calcular los autovalores... ";
            eigVal = A.eigenValues(20, opt);
            eigVal = eigVal.organize();
            if ((((bc &BCL_df) and (bc &BCR_df)) != 0) || (((bc &BCB_df) and (bc &BCT_df)) != 0))
                eigVal(0) = 0.0;
            cout << "Terminado.\n";
            ofstream out(path.data());
            out.precision(15);
            out << sqrt(eigVal)*length/h;
            cout << "\tLos datos se han guardado en: " << path << endl;
        } else {
            in >> eigVal;
            eigVal = pow(eigVal*h/length,2);
            cout << "Terminado.\n";
            cout << "\tLos datos se han importado de: " << path << endl;
        }
    } else {
        cout << "\tCalculando autovalores... ";
        eigVal = A.eigenValues(20, opt);
        eigVal = eigVal.organize();
        if ((((bc &BCL_df) and (bc &BCR_df)) != 0) || (((bc &BCB_df) and (bc &BCT_df)) != 0))
            eigVal(0) = 0.0;
        cout << "Terminado.\n";
    }
    
    cout << "\tCalculando autovectores... ";
    if (((((bc &BCL_df) and (bc &BCR_df)) != 0) || (((bc &BCB_df) and (bc &BCT_df)) != 0)) and eigVal(mode-1) == 0)
        solV.ones();
    else
        solV = A.eigenVector(eigVal(mode-1), opt);
    
    if ((bc &BCL_df || bc &BCB_df) != 0)
        sol.set(0, solV);
    else
        sol.set(1, solV);
    
    cout << "Terminado.\n\n";
    
    return sol;
}

vectorEDP_T EDP::eigenVAL_VEC(vectorEDP_T& x, int mode, unsigned char bc)
{
    return eigenVAL_VEC(x, mode, bc, 0);
}

matrixEDP_T EDP::eigenVAL_VEC(vectorEDP_T &x, vectorEDP_T &y, int modeX, int modeY, unsigned char bc, unsigned char opt)
{
    matrixEDP_T mx = zero<EDP_T>(1,x.dim());
    matrixEDP_T my = zero<EDP_T>(y.dim(),1);
    
    unsigned char bcX = 0, bcY = 0;
    if ((bc &BCL_df) != 0)
        bcX |= BCL_df;
    if ((bc &BCL_f) != 0)
        bcX |= BCL_f;
    if ((bc &BCR_df) != 0)
        bcX |= BCR_df;
    if ((bc &BCR_f) != 0)
        bcX |= BCR_f;
    
    if ((bc &BCB_df) != 0)
        bcY |= BCB_df;
    if ((bc &BCB_f) != 0)
        bcY |= BCB_f;
    if ((bc &BCT_df) != 0)
        bcY |= BCT_df;
    if ((bc &BCT_f) != 0)
        bcY |= BCT_f;
    
    mx.setRowV(0, eigenVAL_VEC(x, modeX, bcX, opt));
    my.setColumnV(0, eigenVAL_VEC(y, modeY, bcY, opt));
    
    return my*mx;
}

matrixEDP_T EDP::eigenVAL_VEC(vectorEDP_T &x, vectorEDP_T &y, int modeX, int modeY, unsigned char bc)
{
    return eigenVAL_VEC(x, y, modeX, modeY, bc, 0);
}



//  -- FUNCIONES DE SALIDA --
void EDP::saveDATA(const string path, const string fileName, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& forSave)
{
    string inPath;
    if ((opt& DESKTOP) != 0) {
        const string home = getenv("HOME");
        inPath = home + "/Desktop/";
    } else if ((opt& DOCUMENTS) != 0) {
        const string home = getenv("HOME");
        inPath = home + "/Documents/";
    } else {
        inPath = path;
    }
    
    if (inPath != "") {
        string solution = inPath + fileName;
        
        ofstream outSolution(solution.data(), ios_base::in | ios_base::out | ios_base::app);
        
        outSolution << forSave;
        
        if (initEDP != true) {
            initEDP = true;
            cout << EDPwarning << "saveDATA(...)] - Los datos se han guardado en: " << solution << endl << endl;
        }
        
    } else
        cout << EDPwarning << "saveDATA(...)] - No se ha especificado ninguna ruta para guardar los datos.\n\n";
}

void EDP::saveDATA(const string path, const string fileName, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& forSave)
{
    saveDATA(path, fileName, 0, x, y, forSave);
}

void EDP::saveDATA(const string fileName, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& forSave)
{
    saveDATA("", fileName, opt, x, y, forSave);
}


void EDP::printMATLAB(const string path, const string fileName, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& forPrint)
{
    int n = forPrint.rows(), m = forPrint.columns();
    EDP_T hx = (EDP_T)(x(m-1)-x(0))/(m-1), hy = (EDP_T)(y(n-1)-y(0))/(n-1);
    
    string inPath;
    if ((opt& DESKTOP) != 0) {
        const string home = getenv("HOME");
        inPath = home + "/Desktop/";
    } else if ((opt& DOCUMENTS) != 0) {
        const string home = getenv("HOME");
        inPath = home + "/Documents/";
    } else {
        inPath = path;
    }
    
    if (inPath != "") {
        const string fileNameData = fileName + ".csv";
        const string solution = inPath + fileNameData;
        
        ofstream outSolution(solution.data());
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < (m-1); j++) {
                outSolution << forPrint(i,j) << ",";
            }
            outSolution << forPrint(i,m-1) << endl;
        }
        
        cout << " Los datos se han guardado en: " << solution << endl << endl;
        
        const string fileNamePlot = "Plot_MATLAB.m";
        const string plotMATLAB = inPath + fileNamePlot;
        
        ofstream outPlot(plotMATLAB.data());
        outPlot << "%%  Plot solutions Laplace/Poisson equation %%\n\n";
        outPlot << "clear all;\n\n";
        outPlot << "%   Importa los datos del archivo\n";
        outPlot << "data = importdata('" << solution << "',',');\n\n";
        outPlot << "%   Crea una maya\n";
        outPlot << "[X,Y] = meshgrid(" << x(0) << ":" << setprecision(4) << hx << ":" << x(m-1) << "," << y(0) << ":" << setprecision(4) << hy << ":" << y(n-1) << ");\n\n";
        outPlot << "%   Dibuja los datos\n";
        outPlot << "subplot(1,2,1)\n";
        outPlot << "meshc(X,Y,data)\n\n";
        outPlot << "xlabel('X = " << x(m-1) << "'); ylabel('Y = " << y(n-1) << "'); zlabel('U(X, Y)');\n";
        outPlot << "title('Laplace/Poisson equation');\n";
        outPlot << "%axis([" << x(0) << x(m-1) << y(0) << " " << y(n-1) << " min(min(data))*1.1 max(max(data))]);\n\n";
        outPlot << "subplot(1,2,2)\n";
        outPlot << "surf(X,Y,data)\n";
        outPlot << "shading interp\n\n";
        outPlot << "xlabel('X = " << x(m-1) << "'); ylabel('Y = " << y(n-1) << "'); zlabel('U(X, Y)');\n";
        outPlot << "title('Laplace/Poisson equation');\n";
        outPlot << "%axis([" << x(0) << x(m-1) << y(0) << " " << y(n-1) << " min(min(data))*1.1 max(max(data))]);\n\n";
        
        cout << " Sin mover de lugar los archivos generados ejecutar el archivo: " << fileNamePlot << " en MATLAB.\n\n";
    } else {
        cout << " [Información printMATLAB]: No se ha especificado ninguna ruta para guardar los datos.\n\n";
    }
}
void EDP::printMATLAB(const string fileName, unsigned char opt, vectorEDP_T& x, vectorEDP_T& y, matrixEDP_T& forPrint)
{
    printMATLAB("", fileName, opt, x, y, forPrint);
}
