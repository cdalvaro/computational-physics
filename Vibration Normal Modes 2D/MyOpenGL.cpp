//
//  MyOpenGL.cpp
//  Vibration Normal Modes 2D
//
//  Created by Carlos David on 5/8/12.
//  Copyright (c) 2012 NEPTUNO. All rights reserved.
//

#include "MyOpenGL.h"
#include "Vectorial.cpp"

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#elif __unix__
#include <GL/gl.h>
#include <GL/glut.h>
#endif

OpenGL_Class myOpenGL;
bool initOpenGL = false;
bool defaultSettings = false;
int plot_number = 1;
int pKeys = 0;
float dAng = 0.5, dPos = 0.1, dZoom = 0.01, zoom = 1.0;
float eyeInitPosX = 1.0, eyeInitPosY = 1.0, eyeInitPosZ = 1.0;
float eyePosX = 1.0, eyePosY = 1.0, eyePosZ = 1.0;

void upPlot()
{
    
}


//  OPENGL - CLASS
//  Init MyOpenGL_Class
void OpenGL_Class::initRendering()
{
    glEnable(GL_DEPTH_TEST);
}


//  Set data
void OpenGL_Class::setData(vector<double> &x, vector<double> &y, unsigned char setDefault)
{
    myOpenGL._vX = x;
    myOpenGL._vY = y;
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Rot_D | Zoom_D | WN_D | dt_D | WSize_D | Update_D);
        defaultSettings = true;
    }
}

void OpenGL_Class::setData(vector<double> &x, vector<double> &y)
{
    setData(x, y, SetDefault);
}

void OpenGL_Class::setData(vector<double> &x, vector<double> &y, vector<double> &z, unsigned char setDefault)
{
    myOpenGL._vX = x;
    myOpenGL._vY = y;
    myOpenGL._vZ = z;
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Rot_D | Zoom_D | WN_D | dt_D | WSize_D | Update_D);
        defaultSettings = true;
    }
}

void OpenGL_Class::setData(vector<double> &x, vector<double> &y, vector<double> &z)
{
    setData(x, y, z, SetDefault);
}

void OpenGL_Class::setData(vector<double> &x, vector<double> &y, matrix<double> &z, unsigned char setDefault)
{
    myOpenGL._vX = x;
    myOpenGL._vY = y;
    myOpenGL._mZ = z;
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Rot_D | Zoom_D | WN_D | dt_D | WSize_D | Update_D);
        defaultSettings = true;
    }
}

void OpenGL_Class::setData(vector<double> &x, vector<double> &y, matrix<double> &z)
{
    setData(x, y, z, SetDefault);
}


//  Set ajustments
void OpenGL_Class::setColors(float red, float green, float blue)
{
    myOpenGL._red = red;
    myOpenGL._green = green;
    myOpenGL._blue = blue;
}

void OpenGL_Class::setZoomsD(float zoomX, float zoomY, unsigned char setDefault)
{
    myOpenGL._zoomX = eyeInitPosX = eyePosX = zoomX;
    myOpenGL._zoomY = eyeInitPosY = eyePosY = zoomY;
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Rot_D | WN_D | dt_D | WSize_D | Update_D);
        defaultSettings = true;
    }
}

void OpenGL_Class::setZooms(float zoomX, float zoomY)
{
    setZoomsD(zoomX, zoomY, SetDefault);
}

void OpenGL_Class::setZoomsD(float zoomX, float zoomY, float zoomZ, unsigned char setDefault)
{
    myOpenGL._zoomX = eyeInitPosX = eyePosX = zoomX;
    myOpenGL._zoomY = eyeInitPosY = eyePosY = zoomY;
    myOpenGL._zoomZ = eyeInitPosZ = eyePosZ = zoomZ;
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Rot_D | WN_D | dt_D | WSize_D | Update_D);
        defaultSettings = true;
    }
}

void OpenGL_Class::setZooms(float zoomX, float zoomY, float zoomZ)
{
    setZoomsD(zoomX, zoomY, zoomZ, SetDefault);
}

void OpenGL_Class::setAxisDirection(float rotX, float rotY, float rotZ, unsigned char setDefault)
{
    myOpenGL._rotX = rotX;
    myOpenGL._rotY = rotY;
    myOpenGL._rotZ = rotZ;
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Zoom_D | WN_D | dt_D | WSize_D | Update_D);
        defaultSettings = true;
    }
}

void OpenGL_Class::setAxisDirection(float rotX, float rotY, float rotZ)
{
    setAxisDirection(rotX, rotY, rotZ, SetDefault);
}

void OpenGL_Class::setStepTime(float dt_milli, unsigned char setDefault)
{
    if (dt_milli >= 1E-4) {
        myOpenGL._dt = (int)ceil(dt_milli*1E4);
        myOpenGL._plots_calc = (int)1;
    } else {
        myOpenGL._dt = (int)1;
        myOpenGL._plots_calc = (int)ceil(1/(dt_milli*1E4));
    }
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Rot_D | Zoom_D | WN_D | WSize_D | Update_D);
        defaultSettings = true;
    }
}

void OpenGL_Class::setStepTime(float dt_milli)
{
    setStepTime(dt_milli, SetDefault);
}

void OpenGL_Class::setWindowSize(float Width, float Height, unsigned char setDefault)
{
    myOpenGL._W = Width;
    myOpenGL._H = Height;
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Rot_D | Zoom_D | WN_D | dt_D | Update_D);
        defaultSettings = true;
    }
}

void OpenGL_Class::setWindowSize(float Width, float Height)
{
    setWindowSize(Width, Height, SetDefault);
}

void OpenGL_Class::setUpdateData(void updateData(), unsigned char setDefault)
{
    myOpenGL._updateData = updateData;
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Rot_D | Zoom_D | WN_D | dt_D | WSize_D);
        defaultSettings = true;
    }
}

void OpenGL_Class::setUpdateData(void updateData())
{
    setUpdateData(updateData, SetDefault);
}

void OpenGL_Class::setWindowName(string windowName, unsigned char setDefault)
{
    myOpenGL._windowName = windowName;
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Rot_D | Zoom_D | dt_D | WSize_D | Update_D);
        defaultSettings = true;
    }
}

void OpenGL_Class::setWindowName(string windowName)
{
    setWindowName(windowName, SetDefault);
}

void OpenGL_Class::setDefaultParameters(unsigned char params)
{
    if ((params& Zoom_D) != 0 && !defaultSettings)
        setZoomsD(1.0, 1.0, 1.0, 0);
    
    if ((params& Rot_D) != 0 && !defaultSettings)
        setAxisDirection(0.0, 0.0, 0.0, 0);
    
    if ((params& WN_D) != 0 && !defaultSettings)
        setWindowName("OpenGL Plot", 0);
    
    if ((params& dt_D) != 0 && !defaultSettings)
        setStepTime(9.1E-4, 0);
    
    if ((params& WSize_D) != 0 && !defaultSettings)
        setWindowSize(800, 800, 0);
    
    if ((params& Update_D) != 0 && !defaultSettings)
        setUpdateData(upPlot, 0);
    
    defaultSettings = true;
}


//  Plot functions
void OpenGL_Class::reshape2D(int Width, int Height)
{
    glViewport (0, 0, (GLsizei)Width, (GLsizei)Height);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    gluPerspective (90, (GLfloat)Width / (GLfloat)Height, -100.0, 10000.0);
    myOpenGL._W = Width;
    myOpenGL._H = Height;
    glMatrixMode (GL_MODELVIEW);
}

void OpenGL_Class::reshape3D(int Width, int Height)
{
    glViewport (0, 0, (GLsizei)Width, (GLsizei)Height);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    gluPerspective(45.0f, Width/Height, 0.1f, 100.0f);
    myOpenGL._W = Width;
    myOpenGL._H = Height;
    glMatrixMode (GL_MODELVIEW);
}

void OpenGL_Class::refreshCameraPosition(int value)
{
    if ((pKeys &UP) != 0)
        myOpenGL._rotX += dAng;
    
    if ((pKeys &DOWN) != 0)
        myOpenGL._rotX -= dAng;
    
    if ((pKeys &LEFT) != 0)
        myOpenGL._rotY -= dAng;
    
    if ((pKeys &RIGHT) != 0)
        myOpenGL._rotY += dAng;
    
    if ((pKeys &Z_in) != 0) {
        eyePosX -= eyePosX*dZoom;
        eyePosY -= eyePosY*dZoom;
        eyePosZ -= eyePosZ*dZoom;
        zoom -= zoom*dZoom;
    }
    
    if ((pKeys &Z_out) != 0) {
        eyePosX += eyePosX*dZoom;
        eyePosY += eyePosY*dZoom;
        eyePosZ += eyePosZ*dZoom;
        zoom += zoom*dZoom;
    }
    
    glutTimerFunc(myOpenGL._plots_calc, refreshCameraPosition, 2);
}


//  2D functions
void OpenGL_Class::display2D()
{
    glClearColor (1.0, 1.0, 1.0, 1.0);
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    if (initOpenGL != true)
        center2D();
    
    glPushMatrix();
	
    myOpenGL.plotAxes2D();
    myOpenGL.drawSolution2D();
    
    glFlush();
    glutSwapBuffers();
    
    initOpenGL = true;
}

void OpenGL_Class::drawSolution2D()
{
    myOpenGL._vMax = myOpenGL._vY.max();
    myOpenGL._vMin = myOpenGL._vY.min();
    
    float adjX, adjY;
    adjX = myOpenGL._cX - myOpenGL._left_right;
    //adjY = myOpenGL._cY - myOpenGL._up_down;
    adjY = - myOpenGL._up_down;
    
    if ((myOpenGL._options& Lines) != 0) {
        glLineWidth(1.0);
        glBegin(GL_LINE_STRIP);
    } else if ((myOpenGL._options& Points) != 0) {
        glBegin(GL_POINTS);
    } else {
        glLineWidth(1.0);
        glBegin(GL_LINE_STRIP);
    }
    
    for (int i=0; i<myOpenGL._vX.dim(); i++) {
        ColorMap::setColor(myOpenGL._vY(i), myOpenGL._vMin, myOpenGL._vMax);
        glColor3f(myOpenGL._red, myOpenGL._green, myOpenGL._blue);
        glVertex2f((myOpenGL._vX(i) - adjX)*myOpenGL._zoomX, (myOpenGL._vY(i) - adjY)*myOpenGL._zoomY);
    }
    
    glEnd();
}

void OpenGL_Class::timerFunction2D(int value)
{
    myOpenGL._updateData();
    
    if (myOpenGL._plots_calc > 1) {
        if (plot_number%myOpenGL._plots_calc == 0){
            plot_number = (int)0;
            myOpenGL.display2D();
        }
    } else
        myOpenGL.display2D();
    
    glutPostRedisplay();
    glutTimerFunc(myOpenGL._dt, timerFunction2D, 0);
}

void OpenGL_Class::plotAxes2D()
{
    const float tmpXmax = myOpenGL._vX.max(), tmpYmax = myOpenGL._vY.max();
    const float tmpLimXmax = myOpenGL._limXmax, tmpLimYmax = myOpenGL._limYmax;
    
    const float tmpXmin = myOpenGL._vX.min(), tmpYmin = myOpenGL._vY.min();
    const float tmpLimXmin = myOpenGL._limXmin, tmpLimYmin = myOpenGL._limYmin;
    
    myOpenGL._limXmax = max(tmpXmax, tmpLimXmax);
    myOpenGL._limYmax = max(tmpYmax, tmpLimYmax);
    
    myOpenGL._limXmin = min(tmpXmin, tmpLimXmin);
    myOpenGL._limYmin = min(tmpYmin, tmpLimYmin);
    
    float compX, compY;
    compX = myOpenGL._cX - myOpenGL._left_right;
    //compY = myOpenGL._cY - myOpenGL._up_down;
    compY = - myOpenGL._up_down;
    
    glLineWidth(1.0);
    glBegin(GL_LINES);
    glColor3f(0.0, 0.0, 0.0);
    glVertex2f((myOpenGL._limXmin*0.01 - compX)*myOpenGL._zoomX, -compY*myOpenGL._zoomY);
    glVertex2f((myOpenGL._limXmax*1.01 - compX)*myOpenGL._zoomX, -compY*myOpenGL._zoomY);
    glVertex2f(-compX*myOpenGL._zoomX, (myOpenGL._limYmin*1.01 - compY)*myOpenGL._zoomY);
    glVertex2f(-compX*myOpenGL._zoomX, (myOpenGL._limYmax*1.01 - compY)*myOpenGL._zoomY);
    glEnd();
}

void OpenGL_Class::center2D()
{
    myOpenGL._cX = (myOpenGL._vX.max() - myOpenGL._vX.min())/2.0;
    myOpenGL._cY = (myOpenGL._vY.max() - myOpenGL._vY.min())/2.0;
}


//  3D functions
void OpenGL_Class::display3D()
{
    glClearColor (1.0f, 1.0f, 1.0f, 1.0f);
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    if (initOpenGL != true)
        myOpenGL.center3D();
    
    gluLookAt(eyePosX, eyePosZ, eyePosY,
              0.0f, 0.0f, 0.0f,
              0.0f, 1.0f, 0.0f);
    
    glRotatef(myOpenGL._rotX, -1.0f, 0.0f, 1.0f);
    glRotatef(myOpenGL._rotY, 0.0f, 1.0f, 0.0f);
    glPushMatrix();
	
    if ((myOpenGL._options &Axis_On) != 0)
        myOpenGL.plotAxes3D();
    
    myOpenGL.drawSolution3D();
    
    glFlush();
    glutSwapBuffers();
    
    initOpenGL = true;
}

void OpenGL_Class::drawSolution3D()
{
    myOpenGL._vMax = myOpenGL._mZ.max();
    myOpenGL._vMin = myOpenGL._mZ.min();
    
    float adjX, adjY, adjZ;
    adjX = -myOpenGL._cX;
    adjY = myOpenGL._left_right - myOpenGL._cY;
    adjZ = myOpenGL._up_down;
    //adjZ = myOpenGL._up_down - myOpenGL._cZ;
    
    if ((myOpenGL._options &Lines) != 0) {
        for (int i=0; i<myOpenGL._vY.dim(); i++) {
            glBegin(GL_LINE_STRIP);
            for (int j=0; j<myOpenGL._vX.dim(); j++) {
                ColorMap::setColor(myOpenGL._mZ(i,j), myOpenGL._vMin, myOpenGL._vMax);
                glColor3f(myOpenGL._red, myOpenGL._green, myOpenGL._blue);
                glVertex3f(myOpenGL._zoomX*(myOpenGL._vX(j) + adjX), myOpenGL._zoomZ*(myOpenGL._mZ(i,j) + adjZ), myOpenGL._zoomY*(myOpenGL._vY(i) + adjY));
            }
            glEnd();
        }
        
        for (int i=0; i<myOpenGL._vX.dim(); i++) {
            glBegin(GL_LINE_STRIP);
            for (int j=0; j<myOpenGL._vY.dim(); j++) {
                ColorMap::setColor(myOpenGL._mZ(j,i), myOpenGL._vMin, myOpenGL._vMax);
                glColor3f(myOpenGL._red, myOpenGL._green, myOpenGL._blue);
                glVertex3f(myOpenGL._zoomX*(myOpenGL._vX(i) + adjX), myOpenGL._zoomZ*(myOpenGL._mZ(j,i) + adjZ), myOpenGL._zoomY*(myOpenGL._vY(j) + adjY));
            }
            glEnd();
        }
    } else if ((myOpenGL._options &Shading) != 0) {
        for (int i=0; i<myOpenGL._vY.dim()-1; i++) {
            glBegin(GL_TRIANGLE_STRIP);
            for (int j=0; j<myOpenGL._vX.dim(); j++) {
                ColorMap::setColor(myOpenGL._mZ(i,j), myOpenGL._vMin, myOpenGL._vMax);
                glColor3f(myOpenGL._red, myOpenGL._green, myOpenGL._blue);
                glVertex3f(myOpenGL._zoomX*(myOpenGL._vX(j) + adjX), myOpenGL._zoomZ*(myOpenGL._mZ(i,j) + adjZ), myOpenGL._zoomY*(myOpenGL._vY(i) + adjY));
                glVertex3f(myOpenGL._zoomX*(myOpenGL._vX(j) + adjX), myOpenGL._zoomZ*(myOpenGL._mZ(i+1,j) + adjZ), myOpenGL._zoomY*(myOpenGL._vY(i+1) + adjY));
            }
            glEnd();
        }
    } else if ((myOpenGL._options &Points) != 0) {
        for (int i=0; i<myOpenGL._vY.dim(); i++) {
            glBegin(GL_POINTS);
            glPointSize(5.0f);
            for (int j=0; j<myOpenGL._vX.dim(); j++) {
                ColorMap::setColor(myOpenGL._mZ(i,j), myOpenGL._vMin, myOpenGL._vMax);
                glColor3f(myOpenGL._red, myOpenGL._green, myOpenGL._blue);
                glVertex3f(myOpenGL._zoomX*(myOpenGL._vX(j) + adjX), myOpenGL._zoomZ*(myOpenGL._mZ(i,j) + adjZ), myOpenGL._zoomY*(myOpenGL._vY(i) + adjY));
            }
            glEnd();
        }
    } else {
        for (int i=0; i<myOpenGL._vY.dim(); i++) {
            glBegin(GL_LINE_STRIP);
            for (int j=0; j<myOpenGL._vX.dim(); j++) {
                ColorMap::setColor(myOpenGL._mZ(i,j), myOpenGL._vMin, myOpenGL._vMax);
                glColor3f(myOpenGL._red, myOpenGL._green, myOpenGL._blue);
                glVertex3f(myOpenGL._zoomX*(myOpenGL._vX(j) + adjX), myOpenGL._zoomZ*(myOpenGL._mZ(i,j) + adjZ), myOpenGL._zoomY*(myOpenGL._vY(i) + adjY));
            }
            glEnd();
        }
        
        for (int i=0; i<myOpenGL._vX.dim(); i++) {
            glBegin(GL_LINE_STRIP);
            for (int j=0; j<myOpenGL._vY.dim(); j++) {
                ColorMap::setColor(myOpenGL._mZ(j,i), myOpenGL._vMin, myOpenGL._vMax);
                glColor3f(myOpenGL._red, myOpenGL._green, myOpenGL._blue);
                glVertex3f(myOpenGL._zoomX*(myOpenGL._vX(i) + adjX), myOpenGL._zoomZ*(myOpenGL._mZ(j,i) + adjZ), myOpenGL._zoomY*(myOpenGL._vY(j) + adjY));
            }
            glEnd();
        }
    }
    
}

void OpenGL_Class::timerFunction3D(int value)
{
    myOpenGL._updateData();
    
    if (myOpenGL._plots_calc > 1) {
        if (plot_number%myOpenGL._plots_calc == 0){
            plot_number = (int)0;
            myOpenGL.display3D();
        }
    } else
        myOpenGL.display3D();
    
    glutPostRedisplay();
    glutTimerFunc(myOpenGL._dt, timerFunction3D, 1);
    plot_number++;
}

void OpenGL_Class::plotAxes3D()
{
    const float tmpXmax = myOpenGL._vX.max(), tmpYmax = myOpenGL._vY.max(), tmpZmax = myOpenGL._mZ.max();
    const float tmpLimXmax = myOpenGL._limXmax, tmpLimYmax = myOpenGL._limYmax, tmpLimZmax = myOpenGL._limZmax;
    
    const float tmpXmin = myOpenGL._vX.min(), tmpYmin = myOpenGL._vY.min(), tmpZmin = myOpenGL._mZ.min();
    const float tmpLimXmin = myOpenGL._limXmin, tmpLimYmin = myOpenGL._limYmin, tmpLimZmin = myOpenGL._limZmin;
    
    myOpenGL._limXmax = max(tmpXmax, tmpLimXmax);
    myOpenGL._limYmax = max(tmpYmax, tmpLimYmax);
    myOpenGL._limZmax = max(tmpZmax, tmpLimZmax);
    
    myOpenGL._limXmin = min(tmpXmin, tmpLimXmin);
    myOpenGL._limYmin = min(tmpYmin, tmpLimYmin);
    myOpenGL._limZmin = min(tmpZmin, tmpLimZmin);
    
    float compX, compY, compZ;
    compX = myOpenGL._cX;
    compY = myOpenGL._cY - myOpenGL._left_right;
    compZ = -myOpenGL._up_down;
    //compZ = myOpenGL._cZ - myOpenGL._up_down;
    
    glLineWidth(1.3);
    glBegin(GL_LINES);
    glColor3f(0.0, 0.0, 0.0);
    glVertex3f(-compX*myOpenGL._zoomX, -compZ*myOpenGL._zoomZ, (myOpenGL._limYmin*1.01-compY)*myOpenGL._zoomY);
    glVertex3f(-compX*myOpenGL._zoomX, -compZ*myOpenGL._zoomZ, (myOpenGL._limYmax*1.01-compY)*myOpenGL._zoomY);
    glVertex3f(-compX*myOpenGL._zoomX, (myOpenGL._limZmin*1.01-compZ)*myOpenGL._zoomZ, -compY*myOpenGL._zoomY);
    glVertex3f(-compX*myOpenGL._zoomX, (myOpenGL._limZmax*1.01-compZ)*myOpenGL._zoomZ, -compY*myOpenGL._zoomY);
    glVertex3f((myOpenGL._limXmin*1.01-compX)*myOpenGL._zoomX, -compZ*myOpenGL._zoomZ, -compY*myOpenGL._zoomY);
    glVertex3f((myOpenGL._limXmax*1.01-compX)*myOpenGL._zoomX, -compZ*myOpenGL._zoomZ, -compY*myOpenGL._zoomY);
    glEnd();
}

void OpenGL_Class::center3D()
{
    if (myOpenGL._vX.max() >= -myOpenGL._vX.min())
        myOpenGL._cX = (myOpenGL._vX.max() - myOpenGL._vX.min())/2.0;
    else
        myOpenGL._cX = -(myOpenGL._vX.max() - myOpenGL._vX.min())/2.0;
    
    if (myOpenGL._vY.max() >= -myOpenGL._vY.min())
        myOpenGL._cY = (myOpenGL._vY.max() - myOpenGL._vY.min())/2.0;
    else
        myOpenGL._cY = -(myOpenGL._vY.max() - myOpenGL._vY.min())/2.0;
    
    if (myOpenGL._mZ.max() >= -myOpenGL._mZ.min())
        myOpenGL._cZ = (myOpenGL._mZ.max() - myOpenGL._mZ.min())/2.0;
    else
        myOpenGL._cZ = -(myOpenGL._mZ.max() - myOpenGL._mZ.min())/2.0;
}


//  Keyboard and mouse functions
void OpenGL_Class::handleKeyPress(unsigned char key, int x, int y)
{
    switch(key)
    {
        case 27:   // escape
            exit(0);
            break;
            
        case '+':
            pKeys |= Z_in;
            break;
            
        case '-':
            pKeys |= Z_out;
            break;
            
        case 'a':
            myOpenGL._left_right -= dPos;
            break;
            
        case 'd':
            myOpenGL._left_right += dPos;
            break;
            
        case 'w':
            myOpenGL._up_down += dPos;
            break;
            
        case 's':
            myOpenGL._up_down -= dPos;
            break;
            
        case 'l':
            myOpenGL._options &= ~Points;
            myOpenGL._options &= ~Shading;
            myOpenGL._options |= Lines;
            break;
            
        case 'k':
            myOpenGL._options &= ~Points;
            myOpenGL._options &= ~Lines;
            myOpenGL._options |= Shading;
            break;
            
        case 'j':
            myOpenGL._options &= ~Lines;
            myOpenGL._options &= ~Shading;
            myOpenGL._options |= Points;
            break;
            
        case 'p':
            eyePosX = (float)(myOpenGL._vX(myOpenGL._vX.dim()-1)-myOpenGL._vX(0))/2;
            eyePosZ = 1.0f;
            eyePosY = (float)(myOpenGL._vY(myOpenGL._vX.dim()-1)-myOpenGL._vX(0))/2;
            myOpenGL._rotX = -35;
            myOpenGL._rotY = -45;
            break;
            
        case '0':
            myOpenGL._left_right = 0.0f;
            myOpenGL._up_down = 0.0f;
            myOpenGL._rotX = 0.0f;
            myOpenGL._rotY = 0.0f;
            myOpenGL._rotZ = 0.0f;
            eyePosX = eyeInitPosX;
            eyePosY = eyeInitPosY;
            eyePosZ = eyeInitPosZ;
            break;
            
        case 'm':
            if ((myOpenGL._options &Axis_On) != 0)
                myOpenGL._options &= ~Axis_On;
            else
                myOpenGL._options |= Axis_On;
            break;
            
        case '1':
            ColorMap::colormap = JET;
            break;
            
        case '2':
            ColorMap::colormap = HOT_to_COLD;
            break;
            
        case '3':
            ColorMap::colormap = COOL;
            break;
            
        case '4':
            ColorMap::colormap = WINTER;
            break;
            
        case '5':
            ColorMap::colormap = HOT;
            break;
            
        case '6':
            ColorMap::colormap = SUMMER;
            break;
            
        case '7':
            ColorMap::colormap = JET_FIXED;
            break;
    }
}

void OpenGL_Class::handleKeyRelease(unsigned char key, int x, int y)
{
    switch (key) {
        case '+':
            pKeys &= ~Z_in;
            break;
        
        case '-':
            pKeys &= ~Z_out;
            break;
    }
}

void OpenGL_Class::handleSpecialKeyPress(int key, int x, int y)
{
    switch (key) {
        case GLUT_KEY_UP:
            pKeys |= UP;
            break;
            
        case GLUT_KEY_DOWN:
            pKeys |= DOWN;
            break;
			
		case GLUT_KEY_LEFT:
            pKeys |= LEFT;
            break;
            
        case GLUT_KEY_RIGHT:
            pKeys |= RIGHT;
            break;
    }
}

void OpenGL_Class::handleSpecialKeyReleased(int key, int x, int y)
{
    switch (key) {
        case GLUT_KEY_UP:
            pKeys &= ~UP;
            break;
            
        case GLUT_KEY_DOWN:
            pKeys &= ~DOWN;
            break;
			
		case GLUT_KEY_LEFT:
            pKeys &= ~LEFT;
            break;
            
        case GLUT_KEY_RIGHT:
            pKeys &= ~RIGHT;
            break;
    }
}


//  COLORMAP - CLASS
float ColorMap_Class::red = 1.0;
float ColorMap_Class::green = 1.0;
float ColorMap_Class::blue = 1.0;
unsigned char ColorMap_Class::colormap = JET;

void ColorMap_Class::Jet(float value, float min, float max)
{
    red = green = blue = 1.0;
    float dv;
    
    if (value < min)
        value = min;
    if (value > max)
        value = max;
    dv = max - min;
    
    if (value < (min + 0.11*dv)) {
        red = 0.0;
        green = 0.0;
        blue = 4.0*(value + 0.8 - min)/dv;
    } else if (value < (min + 0.36*dv)) {
        red = 0.0;
        green = 4.0*(value - min - 0.11*dv)/dv;
    } else if (value < (min + 0.61*dv)) {
        red = 4.0*(value - min - 0.36*dv)/dv;
        blue = 1.0 + 4.0*(min + 0.36*dv - value)/dv;
    } else if (value < (min + 0.86*dv)) {
        blue = 0.0;
        green = 1.0 + 4.0*(min + 0.61*dv - value)/dv;
    } else {
        green = 0.0;
        red = 1.0 + 4.0*(min + 0.86*dv + 0.2 - value)/dv;
        blue = 0.0;
    }
    
    myOpenGL.setColors(red, green, blue);
}

void ColorMap_Class::Hot_to_Cold(float value, float min, float max)
{
    red = green = blue = 1.0;
    float dv;
    
    if (value < min)
        value = min;
    if (value > max)
        value = max;
    dv = max - min;
    
    if (value < (min + 0.25*dv)) {
        red = 0.0;
        green = 4.0*(value - min)/dv;
    } else if (value < (min + 0.5*dv)) {
        red = 0.0;
        blue = 1.0 + 4.0*(min + 0.25*dv - value)/dv;
    } else if (value < (min + 0.75*dv)) {
        red = 4.0*(value - min - 0.5*dv)/dv;
        blue = 0.0;
    } else {
        green = 1.0 + 4.0*(min + 0.75*dv - value)/dv;
        blue = 0.0;
    }
    
    myOpenGL.setColors(red, green, blue);
}

void ColorMap_Class::Cool(float value, float min, float max)
{
    float dv;
    
    if (value < min)
        value = min;
    if (value > max)
        value = max;
    dv = max - min;
    
    red = (value - min)/dv;
    green = 1.0 + (min - value)/dv;
    blue = 1.0;
    
    myOpenGL.setColors(red, green, blue);
}

void ColorMap_Class::Winter(float value, float min, float max)
{
    float dv;
    
    if (value < min)
        value = min;
    if (value > max)
        value = max;
    dv = max - min;
    
    red = 0.0;
    green = (value - min)/dv;
    blue = 1.0 + 0.5*(min - value)/dv;
    
    myOpenGL.setColors(red, green, blue);
}

void ColorMap_Class::Hot(float value, float min, float max)
{
    red = green = blue = 0.0;
    float dv;
    
    if (value < min)
        value = min;
    if (value > max)
        value = max;
    dv = max - min;
    
    if (value < (min + 24.0/64.0*dv)) {
        red = 4.0*(value + 0.2 - min)/dv;
    } else if (value < (min + 48.0/64.0*dv)) {
        red = 1.0;
        green = 4.0*(value - min - 24.0/64.0*dv)/dv;
    } else {
        red = green = 1.0;
        blue = 4.0*(value - min - 48.0/64.0*dv)/dv;
    }
    
    myOpenGL.setColors(red, green, blue);
}

void ColorMap_Class::Summer(float value, float min, float max)
{
    float dv;
    
    if (value < min)
        value = min;
    if (value > max)
        value = max;
    dv = max - min;
    
    red = 1.1*(value - min)/dv;
    green = (value - min + 0.3*dv)/dv;
    blue = 0.4;
    
    myOpenGL.setColors(red, green, blue);
}

void ColorMap_Class::Jet_Fixed(float value, float min, float max)
{
    red = green = blue = 1.0;
    float dv;
    
    if (value < min)
        value = min;
    if (value > max)
        value = max;
    dv = max - min;
    
    if (value < (min + 0.11*dv)) {
        red = 0.0;
        green = 0.0;
        blue = 4.0*(value + 0.8 - min)/dv;
    } else if (value < (min + 0.36*dv)) {
        red = 0.0;
        green = 4.0*(value - min - 0.11*dv)/dv;
    } else if (value < (min + 0.61*dv)) {
        red = 4.0*(value - min - 0.36*dv)/dv;
        blue = 1.0 + 4.0*(min + 0.36*dv - value)/dv;
    } else if (value < (min + 0.86*dv)) {
        blue = 0.0;
        green = 1.0 + 4.0*(min + 0.61*dv - value)/dv;
    } else {
        green = 0.0;
        red = 1.0 + 4.0*(min + 0.86*dv + 0.2 - value)/dv;
        blue = 0.0;
    }
    
    myOpenGL.setColors(red, green, blue);
}

void ColorMap_Class::setColor(float value, float min, float max)
{
    if ((colormap &JET) != 0)
        Jet(value, min, max);
    else if ((colormap &HOT_to_COLD) != 0)
        Hot_to_Cold(value, min, max);
    else if ((colormap &COOL) != 0)
        Cool(value, min, max);
    else if ((colormap &WINTER) != 0)
        Winter(value, min, max);
    else if ((colormap &HOT) != 0)
        Hot(value, min, max);
    else if ((colormap &SUMMER) != 0)
        Summer(value, min, max);
    else if ((colormap &JET_FIXED) != 0)
        Jet_Fixed(value, -0.1, 0.3);
    else
        Jet(value, min, max);
}


//  PLOT FUNCTIONS - NAMESPACE
void OpenGL::plot(vector<double> &x, vector<double> &y, unsigned char options, int argc, const char *argv[])
{
    MyOpenGL_Class::setData(x, y);
    myOpenGL._options = options;
    
    glutInit(&argc, (char **)argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);
    glutInitWindowSize(myOpenGL._W, myOpenGL._H);
    glutInitWindowPosition(100, 100);
    glutCreateWindow((char *)myOpenGL._windowName.c_str());
    MyOpenGL_Class::initRendering();
    glutDisplayFunc(MyOpenGL_Class::display2D);
    glutTimerFunc(myOpenGL._dt, MyOpenGL_Class::timerFunction2D, 0);
    glutTimerFunc(myOpenGL._dt, MyOpenGL_Class::refreshCameraPosition, 0);
    glutKeyboardFunc(MyOpenGL_Class::handleKeyPress);
    glutKeyboardUpFunc(MyOpenGL_Class::handleKeyRelease);
    glutSpecialFunc(MyOpenGL_Class::handleSpecialKeyPress);
    glutSpecialUpFunc(MyOpenGL_Class::handleSpecialKeyReleased);
    glutReshapeFunc(MyOpenGL_Class::reshape2D);
    glutMainLoop();
}

void OpenGL::plot(vector<double> &x, vector<double> &y, vector<double> &z, unsigned char options, int argc, const char *argv[])
{
    MyOpenGL_Class::setData(x, y, z);
    myOpenGL._options = options;
    
    glutInit(&argc, (char **)argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);
    glutInitWindowSize(myOpenGL._W, myOpenGL._H);
    glutInitWindowPosition(100, 100);
    glutCreateWindow((char *)myOpenGL._windowName.c_str());
    MyOpenGL_Class::initRendering();
    glutDisplayFunc(MyOpenGL_Class::display3D);
    glutTimerFunc(myOpenGL._dt, MyOpenGL_Class::timerFunction2D, 0);
    glutTimerFunc(myOpenGL._dt, MyOpenGL_Class::refreshCameraPosition, 0);
    glutKeyboardFunc(MyOpenGL_Class::handleKeyPress);
    glutKeyboardUpFunc(MyOpenGL_Class::handleKeyRelease);
    glutSpecialFunc(MyOpenGL_Class::handleSpecialKeyPress);
    glutSpecialUpFunc(MyOpenGL_Class::handleSpecialKeyReleased);
    glutReshapeFunc(MyOpenGL_Class::reshape3D);
    glutMainLoop();
}

void OpenGL::plot(vector<double> &x, vector<double> &y, matrix<double> &z, unsigned char options, int argc, const char *argv[])
{
    MyOpenGL_Class::setData(x, y, z);
    myOpenGL._options = options;
    
    glutInit(&argc, (char **)argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);
    glutInitWindowSize(myOpenGL._W, myOpenGL._H);
    glutInitWindowPosition(100, 100);
    glutCreateWindow((char *)myOpenGL._windowName.c_str());
    MyOpenGL_Class::initRendering();
    glutDisplayFunc(MyOpenGL_Class::display3D);
    glutTimerFunc(myOpenGL._dt, MyOpenGL_Class::timerFunction3D, 1);
    glutTimerFunc(myOpenGL._dt, MyOpenGL_Class::refreshCameraPosition, 2);
    glutKeyboardFunc(MyOpenGL_Class::handleKeyPress);
    glutKeyboardUpFunc(MyOpenGL_Class::handleKeyRelease);
    glutSpecialFunc(MyOpenGL_Class::handleSpecialKeyPress);
    glutSpecialUpFunc(MyOpenGL_Class::handleSpecialKeyReleased);
    glutReshapeFunc(MyOpenGL_Class::reshape3D);
    glutMainLoop();
}
