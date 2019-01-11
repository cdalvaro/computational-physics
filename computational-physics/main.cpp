//
//  main.cpp
//  Vibration Normal Modes 2D
//
//  Created by Carlos David on 5/8/12.
//  Copyright (c) 2012 NEPTUNO. All rights reserved.
//

#include <iostream>
#include <time.h>

#include "math/differential_equations/SolveEDP.h"
#include "graphics/MyOpenGL.h"

using namespace cda::math::containers;
using namespace cda::math::differential_equations;
using namespace cda::graphics;

namespace cmc = cda::math::containers;

//  Ejemplos
void vibration_modes_A();
void vibration_modes_B();
void vibration_modes_C();
void chladni_patterns_A();
void chladni_patterns_B();
void single_slit_A();
void single_slit_B();
void double_slit();
void cylinder();

#define normalModes     1
#define chladni         2
#define diffraction     3
int model = diffraction;


//  Ejemplos no incluídos en la presentación
void vibration_modes_D();
void pulse();
void chladni_patterns_C();
void single_slit_circular_wave();
void double_slit_circular_wave();
void cylinder_circular_wave();


//  Condiciones de la ecuación de ondas
double QQ(double x, double y);
double BC(double x, double y);


//  Actualización de OpenGL
void calcSol();


//  Función pulso y función sinusoidal
Matrix<double> pick(int x, int y, int rangeX, int rangeY, double strenght);
Matrix<double> sinusoidalForce(int x, int y, int rangeX, int rangeY, double strenght, double freq);


//  Variables compartidas
Matrix<double> cI, cId, sinu;
Matrix<bool> fixedPoints;
Vector<double> vX, vY;
int sPosX, sPosY, sRangeX, sRangeY;
double Ten, P, Lx, Ly, dt, freqSignal, sForce;
double const PI = acos(-1.0);
unsigned char BConditions;
EDP membrane;


int main(int argc, const char * argv[])
{
    std::cout << std::endl;
    std::cout << "   Computación Avanzada\n";
    std::cout << "   Carlos David Álvaro Yunta\n";
    std::cout << "   Proyecto - Ecuación de ondas y cálculo de autovalores/autovectores\n\n";
    
    
    //  Parámetros de la ecuación de ondas: tiempo inicial, valor de c², condiciones de contorno, directorio con autovalores ya calculados
    membrane.time = 0.0;
    membrane.Q2D = QQ;
    membrane.BCT = membrane.BCL = membrane.BCR = membrane.BCB = BC;
    const std::string home = getenv("HOME");
    membrane.pathEDP = home + "/Desktop/Eigen Values/";
    int example;
    
    
    //  ELIGE UNA DE LAS SIGUIENTES TRES OPCIONES
    //  1. normalModes
    //  2. chladni
    //  3. diffraction
    model = normalModes;
    
    //  ELIGE UN EJEMPLO
    //  1. normalModes: 1 - 5
    example = 1;
    
    
    switch (model) {
        case normalModes:
            //  EJEMPLOS - MODOS DE VIBRACIÓN
            switch (example) {
                case 1:
                    vibration_modes_A();            //  Bordes fijos
                    break;
                    
                case 2:
                    vibration_modes_B();            //  Bordes libres
                    break;
                    
                case 3:
                    vibration_modes_C();            //  Composición de modos (bordes fijos)
                    break;
                    
                case 4:
                    vibration_modes_D();            //  Composición de modos (bordes libres)
                    break;
                    
                case 5:
                    pulse();                        //  Pulsos sobre la membrana
                    break;
                    
                default:
                    fprintf(stderr, "Este ejemplo no existe.");
                    break;
            }
            
            //  Plot
            OpenGL::setWindowName("Vibration Modes");
            OpenGL::setStepTime(membrane.dt);
            OpenGL::setUpdateData(calcSol);
            plot(vX, vY, cI, Shading | Axis_On, argc, argv);
            
            break;
            
        case chladni:
            //  EJEMPLOS - FIGURAS DE CHLADNI
            switch (example) {
                case 1:
                    chladni_patterns_A();           //  Una única frecuencia
                    break;
                    
                case 2:
                    chladni_patterns_B();           //  Composición de dos frecuencias
                    break;
                    
                case 3:
                    chladni_patterns_C();           //  Otro ejmemplo de una sola frecuencia
                    break;
                    
                default:
                    fprintf(stderr, "Este ejemplo no existe.");
                    break;
            }
            
            //  Plot
            OpenGL::setWindowName("Chladni Patterns");
            OpenGL::setStepTime(membrane.dt);
            OpenGL::setUpdateData(calcSol);
            plot(vX, vY, cI, Shading, argc, argv);
            break;
            
        case diffraction:
            //  EJEMPLOS - DIFRACCIÓN POR RENDIJAS
            switch (example) {
                case 1:
                    single_slit_A();                //  Una rendija de anchura pequeña
                    break;
                    
                case 2:
                    single_slit_B();                //  Una rendija de anchura grande
                    break;
                    
                case 3:
                    double_slit();                  //  Doble rendija
                    break;
                    
                case 4:
                    cylinder();                     //  Obstáculo: Cilindro
                    break;
                    
                case 5:
                    single_slit_circular_wave();    //  Un rendija con una onda circular
                    break;
                    
                case 6:
                    double_slit_circular_wave();    //  Doble rendija con onda circular
                    break;
                    
                case 7:
                    cylinder_circular_wave();       //  Cilindro con onda circular
                    break;
                    
                default:
                    fprintf(stderr, "Este ejemplo no existe.");
                    break;
            }
            
            //  Plot
            OpenGL::setWindowName("Diffraction");
            OpenGL::setStepTime(membrane.dt);
            OpenGL::setUpdateData(calcSol);
            plot(vX, vY, cI, Shading, argc, argv);
            break;
            
        default:
            fprintf(stderr, "Este modelo no existe.");
            break;
    }
    
    return 0;
}


//  --------------------------------------  //


//  EJEMPLOS - MODOS DE VIBRACIÓN
void vibration_modes_A()                    //  Bordes fijos
{
    //  Difrencial de tiempo
    dt = 0.0005;
    membrane.dt = dt;
    
    //  Propiedades de la membrana
    Ten = 1.0;
    P = 0.1;
    Lx = 1.0;
    Ly = 1.0;
    
    double dx = 0.005, dy = 0.005;
    int n = Ly/dy + 1;
    int m = Lx/dx + 1;
    
    vX = Vector<double>::Zero(m);
    vY = Vector<double>::Zero(n);
    cI = Matrix<double>::Zero(n,m);
    cId = Matrix<double>::Zero(n,m);
    fixedPoints = Matrix<bool>::Zero(n,m);
    
    for (int i=0; i<n; i++)
        vY[i] = i*dy;
    for (int j=0; j<m; j++)
        vX[j] = j*dx;
    
    
    //  CONDICIONES DE CONTORNO
    BConditions = BCT_df | BCB_df | BCL_df | BCR_df;
    
    //  CONDICIÓN INICIAL DE LA VELOCIDAD
    cI = 0.5*membrane.eigenVAL_VEC(vX, vY, 2, 2, BConditions);
}

void vibration_modes_B()                    //  Bordes libres
{
    //  Difrencial de tiempo
    dt = 0.001;
    membrane.dt = dt;
    
    //  Propiedades de la membrana
    Ten = 1.0;
    P = 0.1;
    Lx = 1.0;
    Ly = 1.0;
    
    double dx = 0.005, dy = 0.005;
    int n = Ly/dy + 1;
    int m = Lx/dx + 1;
    
    vX = Vector<double>::Zero(m);
    vY = Vector<double>::Zero(n);
    cI = Matrix<double>::Zero(n,m);
    cId = Matrix<double>::Zero(n,m);
    fixedPoints = Matrix<bool>::Zero(n,m);
    
    for (int i=0; i<n; i++)
        vY[i] = i*dy;
    for (int j=0; j<m; j++)
        vX[j] = j*dx;
    
    
    //  CONDICIONES DE CONTORNO
    BConditions = BCT_df | BCB_df | BCL_df | BCR_df;
    
    //  CONDICIÓN INICIAL DE LA VELOCIDAD
    cI = 0.5*membrane.eigenVAL_VEC(vX, vY, 2, 2, BConditions, SAVE_DATA);
}

void vibration_modes_C()                    //  Composición de modos (bordes fijos)
{
    //  Difrencial de tiempo
    dt = 0.001;
    membrane.dt = dt;
    
    //  Propiedades de la membrana
    Ten = 1.0;
    P = 0.1;
    Lx = 1.0;
    Ly = 1.0;
    
    double dx = 0.005, dy = 0.005;
    int n = Ly/dy + 1;
    int m = Lx/dx + 1;
    
    vX = Vector<double>::Zero(m);
    vY = Vector<double>::Zero(n);
    cI = Matrix<double>::Zero(n,m);
    cId = Matrix<double>::Zero(n,m);
    fixedPoints = Matrix<bool>::Zero(n,m);
    
    for (int i=0; i<n; i++)
        vY[i] = i*dy;
    for (int j=0; j<m; j++)
        vX[j] = j*dx;
    
    
    //  CONDICIONES DE CONTORNO
    BConditions = BCT_f | BCB_f | BCL_f | BCR_f;
    
    //  CONDICIÓN INICIAL DE LA VELOCIDAD
    cI = 0.25*membrane.eigenVAL_VEC(vX, vY, 3, 4, BConditions, IMPORT_DATA);
    cI -= 0.25*membrane.eigenVAL_VEC(vX, vY, 4, 3, BConditions, IMPORT_DATA);
}

void vibration_modes_D()                    //  Composición de modos (bordes libres)
{
    //  Difrencial de tiempo
    dt = 0.001;
    membrane.dt = dt;
    
    //  Propiedades de la membrana
    Ten = 1.0;
    P = 0.1;
    Lx = 1.0;
    Ly = 1.0;
    
    double dx = 0.005, dy = 0.005;
    int n = Ly/dy + 1;
    int m = Lx/dx + 1;
    
    vX = Vector<double>::Zero(m);
    vY = Vector<double>::Zero(n);
    cI = Matrix<double>::Zero(n,m);
    cId = Matrix<double>::Zero(n,m);
    fixedPoints = Matrix<bool>::Zero(n,m);
    
    for (int i=0; i<n; i++)
        vY[i] = i*dy;
    for (int j=0; j<m; j++)
        vX[j] = j*dx;
    
    
    //  CONDICIONES DE CONTORNO
    BConditions = BCT_df | BCB_df | BCL_df | BCR_df;
    
    //  CONDICIÓN INICIAL DE LA VELOCIDAD
    cI = 0.25*membrane.eigenVAL_VEC(vX, vY, 3, 2, BConditions, IMPORT_DATA);
    cI -= 0.25*membrane.eigenVAL_VEC(vX, vY, 2, 3, BConditions, IMPORT_DATA);
}

void pulse()                                //  Pulsos aplicados sobre la membrana
{
    //  Diferencial de tiempo
    dt = 0.0001;
    membrane.dt = dt;
    
    //  Propiedades de la membrana
    Ten = 1.0;
    P = 0.1;
    Lx = 1.0;
    Ly = 1.0;
    
    double dx = 0.005, dy = 0.005;
    int n = Ly/dy + 1;
    int m = Lx/dx + 1;
    
    vX = Vector<double>::Zero(m);
    vY = Vector<double>::Zero(n);
    cI = Matrix<double>::Zero(n,m);
    cId = Matrix<double>::Zero(n,m);
    fixedPoints = Matrix<bool>::Zero(n,m);
    
    for (int i=0; i<n; i++)
        vY[i] = i*dy;
    for (int j=0; j<m; j++)
        vX[j] = j*dx;
    
    
    //  CONDICIONES DE CONTORNO
    BConditions = BCT_df | BCB_df | BCL_df | BCR_df;
    
    //  PULSOS
    cI = pick(m/4, n/4, 3, 3, 0.3);
    cI += pick(3*m/4, 3*n/4, 3, 3, 0.3);
    cI -= pick(m/4, 3*n/4, 3, 3, 0.3);
    cI -= pick(3*m/4, n/4, 3, 3, 0.3);
}


//  EJEMPLOS - FIGURAS DE CHLADNI
void chladni_patterns_A()                   //  Figuras de Chladni
{
    //  Diferencial de tiempo
    dt = 0.0001;
    membrane.dt = dt;
    
    //  Propiedades de la membrana
    Ten = 10.0;
    P = 0.1;
    Lx = 1.0;
    Ly = 1.0;
    
    double dx = 0.005, dy = 0.005;
    int n = Ly/dy + 1;
    int m = Lx/dx + 1;
    
    vX = Vector<double>::Zero(m);
    vY = Vector<double>::Zero(n);
    cI = Matrix<double>::Zero(n,m);
    cId = Matrix<double>::Zero(n,m);
    fixedPoints = Matrix<bool>::Zero(n,m);
    
    for (int i=0; i<n; i++)
        vY[i] = i*dy;
    for (int j=0; j<m; j++)
        vX[j] = j*dx;
    
    
    //  CONDICIONES DE CONTORNO
    BConditions = BCT_df | BCB_df | BCL_df | BCR_df;
    
    //  FRECUENCIA DE ESTIMULACIÓN
    int mx, my;
    mx = 9;             //  Modo X
    my = 7;             //  Modo Y
    
    sPosX = m/2;
    sPosY = n/2;
    sRangeX = 2;
    sRangeY = 2;
    sForce = 0.005;
    freqSignal = (double)sqrt((double)QQ(vX[0],vY[0])*(((double)mx*mx/(Lx*Lx)))+((double)my*my/(Ly*Ly)))*PI;
    
    sinu = Matrix<double>::Zero(n, m);
}

void chladni_patterns_B()                   //  Figuras de Chladni
{
    //  Diferencial de tiempo
    dt = 0.0001;
    membrane.dt = dt;
    
    //  Propiedades de la membrana
    Ten = 10.0;
    P = 0.1;
    Lx = 1.0;
    Ly = 1.0;
    
    double dx = 0.005, dy = 0.005;
    int n = Ly/dy + 1;
    int m = Lx/dx + 1;
    
    vX = Vector<double>::Zero(m);
    vY = Vector<double>::Zero(n);
    cI = Matrix<double>::Zero(n,m);
    cId = Matrix<double>::Zero(n,m);
    fixedPoints = Matrix<bool>::Zero(n,m);
    
    for (int i=0; i<n; i++)
        vY[i] = i*dy;
    for (int j=0; j<m; j++)
        vX[j] = j*dx;
    
    
    //  CONDICIONES DE CONTORNO
    BConditions = BCT_df | BCB_df | BCL_df | BCR_df;
    
    //  FRECUENCIA DE ESTIMULACIÓN
    int mx, my;
    mx = 8;             //  Modo X
    my = 2;             //  Modo Y
    
    sPosX = m/2;
    sPosY = n/2;
    sRangeX = 2;
    sRangeY = 2;
    sForce = 0.005;
    freqSignal = (double)sqrt((double)QQ(vX[0],vY[0])*(((double)mx*mx/(Lx*Lx)))+((double)my*my/(Ly*Ly)))*PI;
    
    mx = 4;             //  Modo X
    my = 6;             //  Modo Y
    freqSignal += (double)sqrt((double)QQ(vX[0],vY[0])*(((double)mx*mx/(Lx*Lx)))+((double)my*my/(Ly*Ly)))*PI;
    
    sinu = Matrix<double>::Zero(n, m);
}

void chladni_patterns_C()                   //  Figuras de Chladni
{
    //  Diferencial de tiempo
    dt = 0.0001;
    membrane.dt = dt;
    
    //  Propiedades de la membrana
    Ten = 10.0;
    P = 0.1;
    Lx = 1.0;
    Ly = 1.0;
    
    double dx = 0.001, dy = 0.001;
    int n = Ly/dy + 1;
    int m = Lx/dx + 1;
    
    vX = Vector<double>::Zero(m);
    vY = Vector<double>::Zero(n);
    cI = Matrix<double>::Zero(n,m);
    cId = Matrix<double>::Zero(n,m);
    fixedPoints = Matrix<bool>::Zero(n,m);
    
    for (int i=0; i<n; i++)
        vY[i] = i*dy;
    for (int j=0; j<m; j++)
        vX[j] = j*dx;
    
    
    //  CONDICIONES DE CONTORNO
    BConditions = BCT_df | BCB_df | BCL_df | BCR_df;
    
    //  FRECUENCIA DE ESTIMULACIÓN
    int mx, my;
    mx = 8;             //  Modo X
    my = 8;             //  Modo Y
    
    sPosX = m/2;
    sPosY = n/2;
    sRangeX = 2;
    sRangeY = 2;
    sForce = 0.005;
    freqSignal = (double)sqrt((double)QQ(vX[0],vY[0])*(((double)mx*mx/(Lx*Lx)))+((double)my*my/(Ly*Ly)))*PI;
    
    sinu = Matrix<double>::Zero(n, m);
}


//  EJEMPLOS - DIFRACCIÓN POR UNA RENDIJA
void single_slit_A()                        //  Una rendija de anchura pequeña
{
    //  Diferencial de tiempo
    dt = 0.0005;
    membrane.dt = dt;
    
    //  Propiedades de la membrana
    Ten = 0.1;
    P = 0.1;
    Lx = 1.0;
    Ly = 1.0;
    
    double dx = 0.005, dy = 0.005;
    int n = Ly/dy + 1;
    int m = Lx/dx + 1;
    
    vX = Vector<double>::Zero(m);
    vY = Vector<double>::Zero(n);
    cI = Matrix<double>::Zero(n,m);
    cId = Matrix<double>::Zero(n,m);
    fixedPoints = Matrix<bool>::Zero(n,m);
    
    for (int i=0; i<n; i++)
        vY[i] = i*dy;
    for (int j=0; j<m; j++)
        vX[j] = j*dx;
    
    
    //  CONDICIONES DE CONTORNO
    BConditions = BCT_df | BCB_df | BCL_df | BCR_df;
    
    //  CONFIGURACIÓN DE LA RENDIJA
    int SlitCenterX = m/3;
    int SlitCenterY = n/2;
    int SlitWidth = n*0.01;
    
    for (int i=0; i<n; i++) {
        fixedPoints[i][SlitCenterX-1] = true;
        fixedPoints[i][SlitCenterX] = true;
        cI[i][SlitCenterX] = 0.5;
        fixedPoints[i][SlitCenterX+1] = true;
    }
    
    for (int i=SlitCenterY-SlitWidth; i<=SlitCenterY+SlitWidth; i++) {
        fixedPoints[i][SlitCenterX] = false;
        fixedPoints[i][SlitCenterX-1] = false;
        fixedPoints[i][SlitCenterX+1] = false;
        cI[i][m/3] = 0.0;
    }
    
    cI[SlitCenterY-(SlitWidth+1)][SlitCenterX] = 0.0;
    cI[SlitCenterY+(SlitWidth+1)][SlitCenterX] = 0.0;
    
    //  FUERZA SINUSOIDAL
    sPosX = 0;                  //  Centro de la fuerza en X
    sPosY = n/2;                //  Centro de la fuerza en Y
    sRangeX = 0;                //  Número de puntos a izquierda y a derecha (2m + 1) puntos en X sobre los que se aplica la fuerza
    sRangeY = n/2;              //  Número de puntos a arriba y abajo (2n + 1) puntos en Y sobre los que se aplica la fuerza
    sForce = 0.3;               //  Amplitud de la fuerza
    freqSignal = 25.0;          //  Frecuencia de oscilación de la fuerza
    sinu = Matrix<double>::Zero(n, m);
    
    //  PARA EL PLOT
    ColorMap::colormap = JET_FIXED;
}

void single_slit_B()                        //  Una rendija de anchura grande
{
    //  Diferencial de tiempo
    dt = 0.0005;
    membrane.dt = dt;
    
    //  Propiedades de la membrana
    Ten = 0.1;
    P = 0.1;
    Lx = 1.0;
    Ly = 1.0;
    
    double dx = 0.005, dy = 0.005;
    int n = Ly/dy + 1;
    int m = Lx/dx + 1;
    
    vX = Vector<double>::Zero(m);
    vY = Vector<double>::Zero(n);
    cI = Matrix<double>::Zero(n,m);
    cId = Matrix<double>::Zero(n,m);
    fixedPoints = Matrix<bool>::Zero(n,m);
    
    for (int i=0; i<n; i++)
        vY[i] = i*dy;
    for (int j=0; j<m; j++)
        vX[j] = j*dx;
    
    
    //  CONDICIONES DE CONTORNO
    BConditions = BCT_df | BCB_df | BCL_df | BCR_df;
    
    //  CONFIGURACIÓN DE LA RENDIJA
    for (int i=0; i<n; i++) {
        fixedPoints[i][m/3-1] = true;
        fixedPoints[i][m/3] = true;
        cI[i][m/3] = 0.5;
        fixedPoints[i][m/3+1] = true;
    }
    
    for (int i=n/2-10; i<=n/2+10; i++) {
        fixedPoints[i][m/3] = false;
        fixedPoints[i][m/3-1] = false;
        fixedPoints[i][m/3+1] = false;
        cI[i][m/3] = 0.0;
    }
    
    cI[n/2-11][m/3] = 0.0;
    cI[n/2+11][m/3] = 0.0;
    
    //  FUERZA SINUSOIDAL
    sPosX = 0;                  //  Centro de la fuerza en X
    sPosY = n/2;                //  Centro de la fuerza en Y
    sRangeX = 0;                //  Número de puntos a izquierda y a derecha (2m + 1) puntos en X sobre los que se aplica la fuerza
    sRangeY = n/2;              //  Número de puntos a arriba y abajo (2n + 1) puntos en Y sobre los que se aplica la fuerza
    sForce = 0.3;               //  Amplitud de la fuerza
    freqSignal = 100.0;         //  Frecuencia de oscilació de la fuerza
    sinu = Matrix<double>::Zero(n, m);
    
    //  PARA EL PLOT
    ColorMap::colormap = JET_FIXED;
}

void double_slit()                          //  Doble rendija
{
    //  Diferencial de tiempo
    dt = 0.0005;
    membrane.dt = dt;
    
    //  Propiedades de la membrana
    Ten = 0.1;
    P = 0.1;
    Lx = 2.0;
    Ly = 2.0;
    
    double dx = 0.006, dy = 0.006;
    int n = Ly/dy + 1;
    int m = Lx/dx + 1;
    
    vX = Vector<double>::Zero(m);
    vY = Vector<double>::Zero(n);
    cI = Matrix<double>::Zero(n,m);
    cId = Matrix<double>::Zero(n,m);
    fixedPoints = Matrix<bool>::Zero(n,m);
    
    for (int i=0; i<n; i++)
        vY[i] = i*dy;
    for (int j=0; j<m; j++)
        vX[j] = j*dx;
    
    
    //  CONDICIONES DE CONTORNO
    BConditions = BCT_df | BCB_df | BCL_df | BCR_df;
    
    //  CONFIGURACIÓN DE LAS RENDIJAS
    for (int i=0; i<n; i++) {
        fixedPoints[i][m/3-1] = true;
        fixedPoints[i][m/3] = true;
        cI[i][m/3] = 0.5;
        fixedPoints[i][m/3+1] = true;
    }
    
    for (int i=2*n/5-3; i<=2*n/5+3; i++) {
        fixedPoints[i][m/3] = false;
        fixedPoints[i][m/3-1] = false;
        fixedPoints[i][m/3+1] = false;
        cI[i][m/3] = 0.0;
    }
    
    for (int i=3*n/5-3; i<=3*n/5+3; i++) {
        fixedPoints[i][m/3] = false;
        fixedPoints[i][m/3-1] = false;
        fixedPoints[i][m/3+1] = false;
        cI[i][m/3] = 0.0;
    }
    
    cI[2*n/5-4][m/3] = 0.0;
    cI[2*n/5+4][m/3] = 0.0;
    cI[3*n/5-4][m/3] = 0.0;
    cI[3*n/5+4][m/3] = 0.0;
    
    //  FUERZA SINUSOIDAL
    sPosX = 0;                  //  Centro de la fuerza en X
    sPosY = n/2;                //  Centro de la fuerza en Y
    sRangeX = 0;                //  Número de puntos a izquierda y a derecha (2m + 1) puntos en X sobre los que se aplica la fuerza
    sRangeY = n/2;              //  Número de puntos a arriba y abajo (2n + 1) puntos en Y sobre los que se aplica la fuerza
    sForce = 0.3;               //  Amplitud de la fuerza
    freqSignal = 100.0;         //  Frecuencia de oscilació de la fuerza
    sinu = Matrix<double>::Zero(n, m);
    
    //  PARA EL PLOT
    ColorMap::colormap = JET_FIXED;
}

void cylinder()                             //  Cilindro
{
    //  Diferencial de tiempo
    dt = 0.0001;
    membrane.dt = dt;
    
    //  Propiedades de la membrana
    Ten = 0.1;
    P = 0.1;
    Lx = 1.0;
    Ly = 1.0;
    
    double dx = 0.005, dy = 0.005;
    int n = Ly/dy + 1;
    int m = Lx/dx + 1;
    
    vX = Vector<double>::Zero(m);
    vY = Vector<double>::Zero(n);
    cI = Matrix<double>::Zero(n,m);
    cId = Matrix<double>::Zero(n,m);
    fixedPoints = Matrix<bool>::Zero(n,m);
    
    for (int i=0; i<n; i++)
        vY[i] = i*dy;
    for (int j=0; j<m; j++)
        vX[j] = j*dx;
    
    
    //  CONDICIONES DE CONTORNO
    BConditions = BCT_df | BCB_df | BCL_df | BCR_df;
    
    //  CONFIGURACIÓN DEL CILINDRO
    int radio = 10;
    
    for (int i=n/2-(radio+1); i<=n/2+(radio+1); i++) {
        for (int j=m/3-(radio+1); j<=m/3+(radio+1); j++) {
            if ((i-n/2)*(i-n/2)+(j-m/3)*(j-m/3) <= (radio+1)*(radio+1)) {
                fixedPoints[i][j] = true;
                cI[i][j] = 0.0;
            }
        }
    }
    
    for (int i=n/2-radio; i<=n/2+radio; i++) {
        for (int j=m/3-radio; j<=m/3+radio; j++) {
            if ((i-n/2)*(i-n/2)+(j-m/3)*(j-m/3) <= radio*radio) {
                fixedPoints[i][j] = true;
                cI[i][j] = 0.5;
            }
        }
    }
    
    //  FUERZA SINUSOIDAL
    sPosX = 0;                  //  Centro de la fuerza en X
    sPosY = n/2;                //  Centro de la fuerza en Y
    sRangeX = 0;                //  Número de puntos a izquierda y a derecha (2m + 1) puntos en X sobre los que se aplica la fuerza
    sRangeY = n/2;              //  Número de puntos a arriba y abajo (2n + 1) puntos en Y sobre los que se aplica la fuerza
    sForce = 0.3;               //  Amplitud de la fuerza
    freqSignal = 100.0;         //  Frecuencia de oscilació de la fuerza
    sinu = Matrix<double>::Zero(n, m);
    
    //  PARA EL PLOT
    ColorMap::colormap = JET_FIXED;
}

void single_slit_circular_wave()            //  Una rendija de anchura con ondas circulares
{
    //  Diferencial de tiempo
    dt = 0.0005;
    membrane.dt = dt;
    
    //  Propiedades de la membrana
    Ten = 0.1;
    P = 0.1;
    Lx = 1.0;
    Ly = 1.0;
    
    double dx = 0.005, dy = 0.005;
    int n = Ly/dy + 1;
    int m = Lx/dx + 1;
    
    vX = Vector<double>::Zero(m);
    vY = Vector<double>::Zero(n);
    cI = Matrix<double>::Zero(n,m);
    cId = Matrix<double>::Zero(n,m);
    fixedPoints = Matrix<bool>::Zero(n,m);
    
    for (int i=0; i<n; i++)
        vY[i] = i*dy;
    for (int j=0; j<m; j++)
        vX[j] = j*dx;
    
    
    //  CONDICIONES DE CONTORNO
    BConditions = BCT_df | BCB_df | BCL_df | BCR_df;
    
    //  CONFIGURACIÓN DE LA RENDIJA
    int SlitCenterX = m/3;
    int SlitCenterY = n/2;
    int SlitWidth = n*0.01;
    
    for (int i=0; i<n; i++) {
        fixedPoints[i][SlitCenterX-1] = true;
        fixedPoints[i][SlitCenterX] = true;
        cI[i][SlitCenterX] = 0.5;
        fixedPoints[i][SlitCenterX+1] = true;
    }
    
    for (int i=SlitCenterY-SlitWidth; i<=SlitCenterY+SlitWidth; i++) {
        fixedPoints[i][SlitCenterX] = false;
        fixedPoints[i][SlitCenterX-1] = false;
        fixedPoints[i][SlitCenterX+1] = false;
        cI[i][m/3] = 0.0;
    }
    
    cI[SlitCenterY-(SlitWidth+1)][SlitCenterX] = 0.0;
    cI[SlitCenterY+(SlitWidth+1)][SlitCenterX] = 0.0;
    
    //  FUERZA SINUSOIDAL
    sPosX = m/6;                //  Centro de la fuerza en X
    sPosY = n/2;                //  Centro de la fuerza en Y
    sRangeX = 2;                //  Número de puntos a izquierda y a derecha (2m + 1) puntos en X sobre los que se aplica la fuerza
    sRangeY = 2;                //  Número de puntos a arriba y abajo (2n + 1) puntos en Y sobre los que se aplica la fuerza
    sForce = 0.3;               //  Amplitud de la fuerza
    freqSignal = 100.0;         //  Frecuencia de oscilació de la fuerza
    sinu = Matrix<double>::Zero(n, m);
    
    //  PARA EL PLOT
    ColorMap::colormap = JET_FIXED;
}

void double_slit_circular_wave()            //  Doble rendija con ondas circulares
{
    //  Diferencial de tiempo
    dt = 0.0005;
    membrane.dt = dt;
    
    //  Propiedades de la membrana
    Ten = 0.1;
    P = 0.1;
    Lx = 2.0;
    Ly = 2.0;
    
    double dx = 0.006, dy = 0.006;
    int n = Ly/dy + 1;
    int m = Lx/dx + 1;
    
    vX = Vector<double>::Zero(m);
    vY = Vector<double>::Zero(n);
    cI = Matrix<double>::Zero(n,m);
    cId = Matrix<double>::Zero(n,m);
    fixedPoints = Matrix<bool>::Zero(n,m);
    
    for (int i=0; i<n; i++)
        vY[i] = i*dy;
    for (int j=0; j<m; j++)
        vX[j] = j*dx;
    
    
    //  CONDICIONES DE CONTORNO
    BConditions = BCT_df | BCB_df | BCL_df | BCR_df;
    
    //  CONFIGURACIÓN DE LAS RENDIJAS
    for (int i=0; i<n; i++) {
        fixedPoints[i][m/3-1] = true;
        fixedPoints[i][m/3] = true;
        cI[i][m/3] = 0.5;
        fixedPoints[i][m/3+1] = true;
    }
    
    for (int i=2*n/5-3; i<=2*n/5+3; i++) {
        fixedPoints[i][m/3] = false;
        fixedPoints[i][m/3-1] = false;
        fixedPoints[i][m/3+1] = false;
        cI[i][m/3] = 0.0;
    }
    
    for (int i=3*n/5-3; i<=3*n/5+3; i++) {
        fixedPoints[i][m/3] = false;
        fixedPoints[i][m/3-1] = false;
        fixedPoints[i][m/3+1] = false;
        cI[i][m/3] = 0.0;
    }
    
    cI[2*n/5-4][m/3] = 0.0;
    cI[2*n/5+4][m/3] = 0.0;
    cI[3*n/5-4][m/3] = 0.0;
    cI[3*n/5+4][m/3] = 0.0;
    
    //  FUERZA SINUSOIDAL
    sPosX = m/6;                //  Centro de la fuerza en X
    sPosY = n/2;                //  Centro de la fuerza en Y
    sRangeX = 2;                //  Número de puntos a izquierda y a derecha (2m + 1) puntos en X sobre los que se aplica la fuerza
    sRangeY = 2;                //  Número de puntos a arriba y abajo (2n + 1) puntos en Y sobre los que se aplica la fuerza
    sForce = 0.3;               //  Amplitud de la fuerza
    freqSignal = 100.0;         //  Frecuencia de oscilació de la fuerza
    sinu = Matrix<double>::Zero(n, m);
    
    //  PARA EL PLOT
    ColorMap::colormap = JET_FIXED;
}

void cylinder_circular_wave()               //  Cilindro con onda circular
{
    //  Diferencial de tiempo
    dt = 0.0001;
    membrane.dt = dt;
    
    //  Propiedades de la membrana
    Ten = 0.1;
    P = 0.1;
    Lx = 1.0;
    Ly = 1.0;
    
    double dx = 0.005, dy = 0.005;
    int n = Ly/dy + 1;
    int m = Lx/dx + 1;
    
    vX = Vector<double>::Zero(m);
    vY = Vector<double>::Zero(n);
    cI = Matrix<double>::Zero(n,m);
    cId = Matrix<double>::Zero(n,m);
    fixedPoints = Matrix<bool>::Zero(n,m);
    
    for (int i=0; i<n; i++)
        vY[i] = i*dy;
    for (int j=0; j<m; j++)
        vX[j] = j*dx;
    
    
    //  CONDICIONES DE CONTORNO
    BConditions = BCT_df | BCB_df | BCL_df | BCR_df;
    
    //  CONFIGURACIÓN DEL CILINDRO
    int radio = 10;
    
    for (int i=n/2-(radio+1); i<=n/2+(radio+1); i++) {
        for (int j=m/2-(radio+1); j<=m/2+(radio+1); j++) {
            if ((i-n/2)*(i-n/2)+(j-m/2)*(j-m/2) <= (radio+1)*(radio+1)) {
                fixedPoints[i][j] = true;
                cI[i][j] = 0.0;
            }
        }
    }
    
    for (int i=n/2-radio; i<=n/2+radio; i++) {
        for (int j=m/2-radio; j<=m/2+radio; j++) {
            if ((i-n/2)*(i-n/2)+(j-m/2)*(j-m/2) <= radio*radio) {
                fixedPoints[i][j] = true;
                cI[i][j] = 0.5;
            }
        }
    }
    
    //  FUERZA SINUSOIDAL
    sPosX = m/4;                //  Centro de la fuerza en X
    sPosY = n/2;                //  Centro de la fuerza en Y
    sRangeX = 2;                //  Número de puntos a izquierda y a derecha (2m + 1) puntos en X sobre los que se aplica la fuerza
    sRangeY = 2;                //  Número de puntos a arriba y abajo (2n + 1) puntos en Y sobre los que se aplica la fuerza
    sForce = 0.3;               //  Amplitud de la fuerza
    freqSignal = 100.0;         //  Frecuencia de oscilació de la fuerza
    sinu = Matrix<double>::Zero(n, m);
    
    //  PARA EL PLOT
    ColorMap::colormap = JET_FIXED;
}


//  --------------------------------------  //


//  CONDICIONES DE LA ECUACIÓN DE ONDAS
double QQ(double x, double y)               //  Valor de c² - [∂²u/∂t² = c²·∆u]
{
    return (Ten * (Lx*Ly) / P);
}

double BC(double x, double y)               //  Condición de contorno
{
    return 0.0;
}


//  ACTUALIZACIÓN DE LA FUNCIÓN
void calcSol()                              //  Calcula los nuevos valores de la membrana y los pinta
{
    if (model == chladni || model == diffraction) {
        cI -= sinu;
        sinu = sinusoidalForce(sPosX, sPosY, sRangeX, sRangeY, sForce, freqSignal);
        cI += sinu;
    }
    
    cI = membrane.solveWave(BConditions, vX, vY, cI, cId, fixedPoints);
    
    if (membrane.dt != dt) {
        OpenGL::setStepTime(membrane.dt);
        dt = membrane.dt;
    }
    
    OpenGL::setData(vX, vY, cI, 0);
}


//  FUNCIÓN PULSO Y FUNCIÓN SINUSOIDAL
Matrix<double> pick(int x, int y, int rangeX, int rangeY, double strenght)                              //  Función pulso
{
    Matrix<double> tmp = Matrix<double>::Zero(cI.Rows(), cI.Columns());
    for (int i=y-rangeY; i<=y+rangeY; i++) {
        for (int j=x-rangeX; j<=x+rangeX; j++) {
            if (i >= 0 && i < cI.Rows() && j >= 0 && j < cI.Columns())
                tmp[i][j] = strenght;
        }
    }
    
    return tmp;
}

Matrix<double> sinusoidalForce(int x, int y, int rangeX, int rangeY, double strenght, double freq)      //  Fuerza sinusoidal
{
    Matrix<double> tmp = Matrix<double>::Zero(cId.Rows(), cId.Columns());
    for (int i=y-rangeY; i<=y+rangeY; i++) {
        for (int j=x-rangeX; j<=x+rangeX; j++) {
            if (i >= 0 && i < cId.Rows() && j >= 0 && j < cId.Columns()) {
                tmp[i][j] = strenght*sin((membrane.time+membrane.dt)*freq);
                fixedPoints[i][j] = true;
            }
        }
    }
    
    return tmp;
}
