//
//  MyOpenGL.cpp
//  Vibration Normal Modes 2D
//
//  Created by Carlos David on 5/8/12.
//  Copyright (c) 2012 NEPTUNO. All rights reserved.
//

#include "MyOpenGL.h"

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#elif __unix__
#include <GL/gl.h>
#include <GL/glut.h>
#endif

using namespace cda::math::containers;
using namespace cda::graphics;

OpenGL OpenGL_Plot;

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
//  Init graphics::Plot
void OpenGL::initRendering()
{
    glEnable(GL_DEPTH_TEST);
}


//  Set data
void OpenGL::setData(const Vector<double> &x, const Vector<double> &y, unsigned char setDefault)
{
    OpenGL_Plot._vX = x;
    OpenGL_Plot._vY = y;
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Rot_D | Zoom_D | WN_D | dt_D | WSize_D | Update_D);
        defaultSettings = true;
    }
}

void OpenGL::setData(const Vector<double> &x, const Vector<double> &y)
{
    setData(x, y, SetDefault);
}

void OpenGL::setData(const Vector<double> &x, const Vector<double> &y, const Vector<double> &z, unsigned char setDefault)
{
    OpenGL_Plot._vX = x;
    OpenGL_Plot._vY = y;
    OpenGL_Plot._vZ = z;
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Rot_D | Zoom_D | WN_D | dt_D | WSize_D | Update_D);
        defaultSettings = true;
    }
}

void OpenGL::setData(const Vector<double> &x, const Vector<double> &y, const Vector<double> &z)
{
    setData(x, y, z, SetDefault);
}

void OpenGL::setData(const Vector<double> &x, const Vector<double> &y, const Matrix<double> &z, unsigned char setDefault)
{
    OpenGL_Plot._vX = x;
    OpenGL_Plot._vY = y;
    OpenGL_Plot._mZ = z;
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Rot_D | Zoom_D | WN_D | dt_D | WSize_D | Update_D);
        defaultSettings = true;
    }
}

void OpenGL::setData(const Vector<double> &x, const Vector<double> &y, const Matrix<double> &z)
{
    setData(x, y, z, SetDefault);
}


//  Set ajustments
void OpenGL::setColors(float red, float green, float blue)
{
    OpenGL_Plot._red = red;
    OpenGL_Plot._green = green;
    OpenGL_Plot._blue = blue;
}

void OpenGL::setZoomsD(float zoomX, float zoomY, unsigned char setDefault)
{
    OpenGL_Plot._zoomX = eyeInitPosX = eyePosX = zoomX;
    OpenGL_Plot._zoomY = eyeInitPosY = eyePosY = zoomY;
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Rot_D | WN_D | dt_D | WSize_D | Update_D);
        defaultSettings = true;
    }
}

void OpenGL::setZooms(float zoomX, float zoomY)
{
    setZoomsD(zoomX, zoomY, SetDefault);
}

void OpenGL::setZoomsD(float zoomX, float zoomY, float zoomZ, unsigned char setDefault)
{
    OpenGL_Plot._zoomX = eyeInitPosX = eyePosX = zoomX;
    OpenGL_Plot._zoomY = eyeInitPosY = eyePosY = zoomY;
    OpenGL_Plot._zoomZ = eyeInitPosZ = eyePosZ = zoomZ;
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Rot_D | WN_D | dt_D | WSize_D | Update_D);
        defaultSettings = true;
    }
}

void OpenGL::setZooms(float zoomX, float zoomY, float zoomZ)
{
    setZoomsD(zoomX, zoomY, zoomZ, SetDefault);
}

void OpenGL::setAxisDirection(float rotX, float rotY, float rotZ, unsigned char setDefault)
{
    OpenGL_Plot._rotX = rotX;
    OpenGL_Plot._rotY = rotY;
    OpenGL_Plot._rotZ = rotZ;
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Zoom_D | WN_D | dt_D | WSize_D | Update_D);
        defaultSettings = true;
    }
}

void OpenGL::setAxisDirection(float rotX, float rotY, float rotZ)
{
    setAxisDirection(rotX, rotY, rotZ, SetDefault);
}

void OpenGL::setStepTime(float dt_milli, unsigned char setDefault)
{
    if (dt_milli >= 1E-4) {
        OpenGL_Plot._dt = (int)ceil(dt_milli*1E4);
        OpenGL_Plot._plots_calc = (int)1;
    } else {
        OpenGL_Plot._dt = (int)1;
        OpenGL_Plot._plots_calc = (int)ceil(1/(dt_milli*1E4));
    }
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Rot_D | Zoom_D | WN_D | WSize_D | Update_D);
        defaultSettings = true;
    }
}

void OpenGL::setStepTime(float dt_milli)
{
    setStepTime(dt_milli, SetDefault);
}

void OpenGL::setWindowSize(float Width, float Height, unsigned char setDefault)
{
    OpenGL_Plot._W = Width;
    OpenGL_Plot._H = Height;
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Rot_D | Zoom_D | WN_D | dt_D | Update_D);
        defaultSettings = true;
    }
}

void OpenGL::setWindowSize(float Width, float Height)
{
    setWindowSize(Width, Height, SetDefault);
}

void OpenGL::setUpdateData(void updateData(), unsigned char setDefault)
{
    OpenGL_Plot._updateData = updateData;
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Rot_D | Zoom_D | WN_D | dt_D | WSize_D);
        defaultSettings = true;
    }
}

void OpenGL::setUpdateData(void updateData())
{
    setUpdateData(updateData, SetDefault);
}

void OpenGL::setWindowName(std::string windowName, unsigned char setDefault)
{
    OpenGL_Plot._windowName = windowName;
    
    if (!defaultSettings && (setDefault& SetDefault) != 0) {
        setDefaultParameters(Rot_D | Zoom_D | dt_D | WSize_D | Update_D);
        defaultSettings = true;
    }
}

void OpenGL::setWindowName(std::string windowName)
{
    setWindowName(windowName, SetDefault);
}

void OpenGL::setDefaultParameters(unsigned char params)
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
void OpenGL::reshape2D(int Width, int Height)
{
    glViewport (0, 0, (GLsizei)Width, (GLsizei)Height);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    gluPerspective (90, (GLfloat)Width / (GLfloat)Height, -100.0, 10000.0);
    OpenGL_Plot._W = Width;
    OpenGL_Plot._H = Height;
    glMatrixMode (GL_MODELVIEW);
}

void OpenGL::reshape3D(int Width, int Height)
{
    glViewport (0, 0, (GLsizei)Width, (GLsizei)Height);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    gluPerspective(45.0f, Width/Height, 0.1f, 100.0f);
    OpenGL_Plot._W = Width;
    OpenGL_Plot._H = Height;
    glMatrixMode (GL_MODELVIEW);
}

void OpenGL::refreshCameraPosition(int value)
{
    if ((pKeys &UP) != 0)
        OpenGL_Plot._rotX += dAng;
    
    if ((pKeys &DOWN) != 0)
        OpenGL_Plot._rotX -= dAng;
    
    if ((pKeys &LEFT) != 0)
        OpenGL_Plot._rotY -= dAng;
    
    if ((pKeys &RIGHT) != 0)
        OpenGL_Plot._rotY += dAng;
    
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
    
    glutTimerFunc(OpenGL_Plot._plots_calc, refreshCameraPosition, 2);
}


//  2D functions
void OpenGL::display2D()
{
    glClearColor (1.0, 1.0, 1.0, 1.0);
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    if (initOpenGL != true)
        center2D();
    
    glPushMatrix();
	
    OpenGL_Plot.plotAxes2D();
    OpenGL_Plot.drawSolution2D();
    
    glFlush();
    glutSwapBuffers();
    
    initOpenGL = true;
}

void OpenGL::drawSolution2D()
{
    OpenGL_Plot._vMax = OpenGL_Plot._vY.MaximumElement();
    OpenGL_Plot._vMin = OpenGL_Plot._vY.MinimumElement();
    
    float adjX, adjY;
    adjX = OpenGL_Plot._cX - OpenGL_Plot._left_right;
    //adjY = OpenGL_Plot._cY - OpenGL_Plot._up_down;
    adjY = - OpenGL_Plot._up_down;
    
    if ((OpenGL_Plot._options& Lines) != 0) {
        glLineWidth(1.0);
        glBegin(GL_LINE_STRIP);
    } else if ((OpenGL_Plot._options& Points) != 0) {
        glBegin(GL_POINTS);
    } else {
        glLineWidth(1.0);
        glBegin(GL_LINE_STRIP);
    }
    
    for (int i=0; i<OpenGL_Plot._vX.Size(); i++) {
        ColorMap::setColor(OpenGL_Plot._vY[i], OpenGL_Plot._vMin, OpenGL_Plot._vMax);
        glColor3f(OpenGL_Plot._red, OpenGL_Plot._green, OpenGL_Plot._blue);
        glVertex2f((OpenGL_Plot._vX[i] - adjX)*OpenGL_Plot._zoomX, (OpenGL_Plot._vY[i] - adjY)*OpenGL_Plot._zoomY);
    }
    
    glEnd();
}

void OpenGL::timerFunction2D(int value)
{
    OpenGL_Plot._updateData();
    
    if (OpenGL_Plot._plots_calc > 1) {
        if (plot_number%OpenGL_Plot._plots_calc == 0){
            plot_number = (int)0;
            OpenGL_Plot.display2D();
        }
    } else
        OpenGL_Plot.display2D();
    
    glutPostRedisplay();
    glutTimerFunc(OpenGL_Plot._dt, timerFunction2D, 0);
}

void OpenGL::plotAxes2D()
{
    const float tmpXmax = OpenGL_Plot._vX.MaximumElement(), tmpYmax = OpenGL_Plot._vY.MaximumElement();
    const float tmpLimXmax = OpenGL_Plot._limXmax, tmpLimYmax = OpenGL_Plot._limYmax;
    
    const float tmpXmin = OpenGL_Plot._vX.MinimumElement(), tmpYmin = OpenGL_Plot._vY.MinimumElement();
    const float tmpLimXmin = OpenGL_Plot._limXmin, tmpLimYmin = OpenGL_Plot._limYmin;
    
    OpenGL_Plot._limXmax = std::max(tmpXmax, tmpLimXmax);
    OpenGL_Plot._limYmax = std::max(tmpYmax, tmpLimYmax);
    
    OpenGL_Plot._limXmin = std::min(tmpXmin, tmpLimXmin);
    OpenGL_Plot._limYmin = std::min(tmpYmin, tmpLimYmin);
    
    float compX, compY;
    compX = OpenGL_Plot._cX - OpenGL_Plot._left_right;
    //compY = OpenGL_Plot._cY - OpenGL_Plot._up_down;
    compY = - OpenGL_Plot._up_down;
    
    glLineWidth(1.0);
    glBegin(GL_LINES);
    glColor3f(0.0, 0.0, 0.0);
    glVertex2f((OpenGL_Plot._limXmin*0.01 - compX)*OpenGL_Plot._zoomX, -compY*OpenGL_Plot._zoomY);
    glVertex2f((OpenGL_Plot._limXmax*1.01 - compX)*OpenGL_Plot._zoomX, -compY*OpenGL_Plot._zoomY);
    glVertex2f(-compX*OpenGL_Plot._zoomX, (OpenGL_Plot._limYmin*1.01 - compY)*OpenGL_Plot._zoomY);
    glVertex2f(-compX*OpenGL_Plot._zoomX, (OpenGL_Plot._limYmax*1.01 - compY)*OpenGL_Plot._zoomY);
    glEnd();
}

void OpenGL::center2D()
{
    OpenGL_Plot._cX = (OpenGL_Plot._vX.MaximumElement() - OpenGL_Plot._vX.MinimumElement())/2.0;
    OpenGL_Plot._cY = (OpenGL_Plot._vY.MaximumElement() - OpenGL_Plot._vY.MinimumElement())/2.0;
}


//  3D functions
void OpenGL::display3D()
{
    glClearColor (1.0f, 1.0f, 1.0f, 1.0f);
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    if (initOpenGL != true)
        OpenGL_Plot.center3D();
    
    gluLookAt(eyePosX, eyePosZ, eyePosY,
              0.0f, 0.0f, 0.0f,
              0.0f, 1.0f, 0.0f);
    
    glRotatef(OpenGL_Plot._rotX, -1.0f, 0.0f, 1.0f);
    glRotatef(OpenGL_Plot._rotY, 0.0f, 1.0f, 0.0f);
    glPushMatrix();
	
    if ((OpenGL_Plot._options &Axis_On) != 0)
        OpenGL_Plot.plotAxes3D();
    
    OpenGL_Plot.drawSolution3D();
    
    glFlush();
    glutSwapBuffers();
    
    initOpenGL = true;
}

void OpenGL::drawSolution3D()
{
    OpenGL_Plot._vMax = OpenGL_Plot._mZ.MaximumElement();
    OpenGL_Plot._vMin = OpenGL_Plot._mZ.MinimumElement();
    
    float adjX, adjY, adjZ;
    adjX = -OpenGL_Plot._cX;
    adjY = OpenGL_Plot._left_right - OpenGL_Plot._cY;
    adjZ = OpenGL_Plot._up_down;
    //adjZ = OpenGL_Plot._up_down - OpenGL_Plot._cZ;
    
    if ((OpenGL_Plot._options &Lines) != 0) {
        for (int i=0; i<OpenGL_Plot._vY.Size(); i++) {
            glBegin(GL_LINE_STRIP);
            for (int j=0; j<OpenGL_Plot._vX.Size(); j++) {
                ColorMap::setColor(OpenGL_Plot._mZ[i][j], OpenGL_Plot._vMin, OpenGL_Plot._vMax);
                glColor3f(OpenGL_Plot._red, OpenGL_Plot._green, OpenGL_Plot._blue);
                glVertex3f(OpenGL_Plot._zoomX*(OpenGL_Plot._vX[j] + adjX), OpenGL_Plot._zoomZ*(OpenGL_Plot._mZ[i][j] + adjZ), OpenGL_Plot._zoomY*(OpenGL_Plot._vY[i] + adjY));
            }
            glEnd();
        }
        
        for (int i=0; i<OpenGL_Plot._vX.Size(); i++) {
            glBegin(GL_LINE_STRIP);
            for (int j=0; j<OpenGL_Plot._vY.Size(); j++) {
                ColorMap::setColor(OpenGL_Plot._mZ[j][i], OpenGL_Plot._vMin, OpenGL_Plot._vMax);
                glColor3f(OpenGL_Plot._red, OpenGL_Plot._green, OpenGL_Plot._blue);
                glVertex3f(OpenGL_Plot._zoomX*(OpenGL_Plot._vX[i] + adjX), OpenGL_Plot._zoomZ*(OpenGL_Plot._mZ[j][i] + adjZ), OpenGL_Plot._zoomY*(OpenGL_Plot._vY[j] + adjY));
            }
            glEnd();
        }
    } else if ((OpenGL_Plot._options &Shading) != 0) {
        for (int i=0; i<OpenGL_Plot._vY.Size()-1; i++) {
            glBegin(GL_TRIANGLE_STRIP);
            for (int j=0; j<OpenGL_Plot._vX.Size(); j++) {
                ColorMap::setColor(OpenGL_Plot._mZ[i][j], OpenGL_Plot._vMin, OpenGL_Plot._vMax);
                glColor3f(OpenGL_Plot._red, OpenGL_Plot._green, OpenGL_Plot._blue);
                glVertex3f(OpenGL_Plot._zoomX*(OpenGL_Plot._vX[j] + adjX), OpenGL_Plot._zoomZ*(OpenGL_Plot._mZ[i][j] + adjZ), OpenGL_Plot._zoomY*(OpenGL_Plot._vY[i] + adjY));
                glVertex3f(OpenGL_Plot._zoomX*(OpenGL_Plot._vX[j] + adjX), OpenGL_Plot._zoomZ*(OpenGL_Plot._mZ[i+1][j] + adjZ), OpenGL_Plot._zoomY*(OpenGL_Plot._vY[i+1] + adjY));
            }
            glEnd();
        }
    } else if ((OpenGL_Plot._options &Points) != 0) {
        for (int i=0; i<OpenGL_Plot._vY.Size(); i++) {
            glBegin(GL_POINTS);
            glPointSize(5.0f);
            for (int j=0; j<OpenGL_Plot._vX.Size(); j++) {
                ColorMap::setColor(OpenGL_Plot._mZ[i][j], OpenGL_Plot._vMin, OpenGL_Plot._vMax);
                glColor3f(OpenGL_Plot._red, OpenGL_Plot._green, OpenGL_Plot._blue);
                glVertex3f(OpenGL_Plot._zoomX*(OpenGL_Plot._vX[j] + adjX), OpenGL_Plot._zoomZ*(OpenGL_Plot._mZ[i][j] + adjZ), OpenGL_Plot._zoomY*(OpenGL_Plot._vY[i] + adjY));
            }
            glEnd();
        }
    } else {
        for (int i=0; i<OpenGL_Plot._vY.Size(); i++) {
            glBegin(GL_LINE_STRIP);
            for (int j=0; j<OpenGL_Plot._vX.Size(); j++) {
                ColorMap::setColor(OpenGL_Plot._mZ[i][j], OpenGL_Plot._vMin, OpenGL_Plot._vMax);
                glColor3f(OpenGL_Plot._red, OpenGL_Plot._green, OpenGL_Plot._blue);
                glVertex3f(OpenGL_Plot._zoomX*(OpenGL_Plot._vX[j] + adjX), OpenGL_Plot._zoomZ*(OpenGL_Plot._mZ[i][j] + adjZ), OpenGL_Plot._zoomY*(OpenGL_Plot._vY[i] + adjY));
            }
            glEnd();
        }
        
        for (int i=0; i<OpenGL_Plot._vX.Size(); i++) {
            glBegin(GL_LINE_STRIP);
            for (int j=0; j<OpenGL_Plot._vY.Size(); j++) {
                ColorMap::setColor(OpenGL_Plot._mZ[j][i], OpenGL_Plot._vMin, OpenGL_Plot._vMax);
                glColor3f(OpenGL_Plot._red, OpenGL_Plot._green, OpenGL_Plot._blue);
                glVertex3f(OpenGL_Plot._zoomX*(OpenGL_Plot._vX[i] + adjX), OpenGL_Plot._zoomZ*(OpenGL_Plot._mZ[j][i] + adjZ), OpenGL_Plot._zoomY*(OpenGL_Plot._vY[j] + adjY));
            }
            glEnd();
        }
    }
    
}

void OpenGL::timerFunction3D(int value)
{
    OpenGL_Plot._updateData();
    
    if (OpenGL_Plot._plots_calc > 1) {
        if (plot_number%OpenGL_Plot._plots_calc == 0){
            plot_number = (int)0;
            OpenGL_Plot.display3D();
        }
    } else
        OpenGL_Plot.display3D();
    
    glutPostRedisplay();
    glutTimerFunc(OpenGL_Plot._dt, timerFunction3D, 1);
    plot_number++;
}

void OpenGL::plotAxes3D()
{
    const float tmpXmax = OpenGL_Plot._vX.MaximumElement(), tmpYmax = OpenGL_Plot._vY.MaximumElement(), tmpZmax = OpenGL_Plot._mZ.MaximumElement();
    const float tmpLimXmax = OpenGL_Plot._limXmax, tmpLimYmax = OpenGL_Plot._limYmax, tmpLimZmax = OpenGL_Plot._limZmax;
    
    const float tmpXmin = OpenGL_Plot._vX.MinimumElement(), tmpYmin = OpenGL_Plot._vY.MinimumElement(), tmpZmin = OpenGL_Plot._mZ.MinimumElement();
    const float tmpLimXmin = OpenGL_Plot._limXmin, tmpLimYmin = OpenGL_Plot._limYmin, tmpLimZmin = OpenGL_Plot._limZmin;
    
    OpenGL_Plot._limXmax = std::max(tmpXmax, tmpLimXmax);
    OpenGL_Plot._limYmax = std::max(tmpYmax, tmpLimYmax);
    OpenGL_Plot._limZmax = std::max(tmpZmax, tmpLimZmax);
    
    OpenGL_Plot._limXmin = std::min(tmpXmin, tmpLimXmin);
    OpenGL_Plot._limYmin = std::min(tmpYmin, tmpLimYmin);
    OpenGL_Plot._limZmin = std::min(tmpZmin, tmpLimZmin);
    
    float compX, compY, compZ;
    compX = OpenGL_Plot._cX;
    compY = OpenGL_Plot._cY - OpenGL_Plot._left_right;
    compZ = -OpenGL_Plot._up_down;
    //compZ = OpenGL_Plot._cZ - OpenGL_Plot._up_down;
    
    glLineWidth(1.3);
    glBegin(GL_LINES);
    glColor3f(0.0, 0.0, 0.0);
    glVertex3f(-compX*OpenGL_Plot._zoomX, -compZ*OpenGL_Plot._zoomZ, (OpenGL_Plot._limYmin*1.01-compY)*OpenGL_Plot._zoomY);
    glVertex3f(-compX*OpenGL_Plot._zoomX, -compZ*OpenGL_Plot._zoomZ, (OpenGL_Plot._limYmax*1.01-compY)*OpenGL_Plot._zoomY);
    glVertex3f(-compX*OpenGL_Plot._zoomX, (OpenGL_Plot._limZmin*1.01-compZ)*OpenGL_Plot._zoomZ, -compY*OpenGL_Plot._zoomY);
    glVertex3f(-compX*OpenGL_Plot._zoomX, (OpenGL_Plot._limZmax*1.01-compZ)*OpenGL_Plot._zoomZ, -compY*OpenGL_Plot._zoomY);
    glVertex3f((OpenGL_Plot._limXmin*1.01-compX)*OpenGL_Plot._zoomX, -compZ*OpenGL_Plot._zoomZ, -compY*OpenGL_Plot._zoomY);
    glVertex3f((OpenGL_Plot._limXmax*1.01-compX)*OpenGL_Plot._zoomX, -compZ*OpenGL_Plot._zoomZ, -compY*OpenGL_Plot._zoomY);
    glEnd();
}

void OpenGL::center3D() {
    
    auto max_vX = OpenGL_Plot._vX.MaximumElement();
    auto min_vX = OpenGL_Plot._vX.MinimumElement();
    OpenGL_Plot._cX = (max_vX >= -min_vX ? 0.5 : -0.5) * (max_vX - min_vX);
    
    auto max_vY = OpenGL_Plot._vY.MaximumElement();
    auto min_vY = OpenGL_Plot._vY.MinimumElement();
    OpenGL_Plot._cY = (max_vY >= -min_vY ? 0.5 : -0.5) * (max_vY - min_vY);
    
    auto max_mZ = OpenGL_Plot._mZ.MaximumElement();
    auto min_mZ = OpenGL_Plot._mZ.MinimumElement();
    OpenGL_Plot._cZ = (max_mZ >= -min_mZ ? 0.5 : -0.5) * (max_mZ - min_mZ);
}


//  Keyboard and mouse functions
void OpenGL::handleKeyPress(unsigned char key, int x, int y)
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
            OpenGL_Plot._left_right -= dPos;
            break;
            
        case 'd':
            OpenGL_Plot._left_right += dPos;
            break;
            
        case 'w':
            OpenGL_Plot._up_down += dPos;
            break;
            
        case 's':
            OpenGL_Plot._up_down -= dPos;
            break;
            
        case 'l':
            OpenGL_Plot._options &= ~Points;
            OpenGL_Plot._options &= ~Shading;
            OpenGL_Plot._options |= Lines;
            break;
            
        case 'k':
            OpenGL_Plot._options &= ~Points;
            OpenGL_Plot._options &= ~Lines;
            OpenGL_Plot._options |= Shading;
            break;
            
        case 'j':
            OpenGL_Plot._options &= ~Lines;
            OpenGL_Plot._options &= ~Shading;
            OpenGL_Plot._options |= Points;
            break;
            
        case 'p':
            eyePosX = (float)(OpenGL_Plot._vX[OpenGL_Plot._vX.Size()-1]-OpenGL_Plot._vX[0])/2;
            eyePosZ = 1.0f;
            eyePosY = (float)(OpenGL_Plot._vY[OpenGL_Plot._vX.Size()-1]-OpenGL_Plot._vX[0])/2;
            OpenGL_Plot._rotX = -35;
            OpenGL_Plot._rotY = -45;
            break;
            
        case '0':
            OpenGL_Plot._left_right = 0.0f;
            OpenGL_Plot._up_down = 0.0f;
            OpenGL_Plot._rotX = 0.0f;
            OpenGL_Plot._rotY = 0.0f;
            OpenGL_Plot._rotZ = 0.0f;
            eyePosX = eyeInitPosX;
            eyePosY = eyeInitPosY;
            eyePosZ = eyeInitPosZ;
            break;
            
        case 'm':
            if ((OpenGL_Plot._options &Axis_On) != 0)
                OpenGL_Plot._options &= ~Axis_On;
            else
                OpenGL_Plot._options |= Axis_On;
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

void OpenGL::handleKeyRelease(unsigned char key, int x, int y)
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

void OpenGL::handleSpecialKeyPress(int key, int x, int y)
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

void OpenGL::handleSpecialKeyReleased(int key, int x, int y)
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
float ColorMap::red = 1.0;
float ColorMap::green = 1.0;
float ColorMap::blue = 1.0;
unsigned char ColorMap::colormap = JET;

void ColorMap::Jet(float value, float min, float max)
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
    
    OpenGL_Plot.setColors(red, green, blue);
}

void ColorMap::Hot_to_Cold(float value, float min, float max)
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
    
    OpenGL_Plot.setColors(red, green, blue);
}

void ColorMap::Cool(float value, float min, float max)
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
    
    OpenGL_Plot.setColors(red, green, blue);
}

void ColorMap::Winter(float value, float min, float max)
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
    
    OpenGL_Plot.setColors(red, green, blue);
}

void ColorMap::Hot(float value, float min, float max)
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
    
    OpenGL_Plot.setColors(red, green, blue);
}

void ColorMap::Summer(float value, float min, float max)
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
    
    OpenGL_Plot.setColors(red, green, blue);
}

void ColorMap::Jet_Fixed(float value, float min, float max)
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
    
    OpenGL_Plot.setColors(red, green, blue);
}

void ColorMap::setColor(float value, float min, float max)
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
void cda::graphics::plot(Vector<double> &x, Vector<double> &y, unsigned char options, int argc, const char *argv[])
{
    OpenGL::setData(x, y);
    OpenGL_Plot._options = options;
    
    glutInit(&argc, (char **)argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);
    glutInitWindowSize(OpenGL_Plot._W, OpenGL_Plot._H);
    glutInitWindowPosition(100, 100);
    glutCreateWindow((char *)OpenGL_Plot._windowName.c_str());
    OpenGL::initRendering();
    glutDisplayFunc(OpenGL::display2D);
    glutTimerFunc(OpenGL_Plot._dt, OpenGL::timerFunction2D, 0);
    glutTimerFunc(OpenGL_Plot._dt, OpenGL::refreshCameraPosition, 0);
    glutKeyboardFunc(OpenGL::handleKeyPress);
    glutKeyboardUpFunc(OpenGL::handleKeyRelease);
    glutSpecialFunc(OpenGL::handleSpecialKeyPress);
    glutSpecialUpFunc(OpenGL::handleSpecialKeyReleased);
    glutReshapeFunc(OpenGL::reshape2D);
    glutMainLoop();
}

void cda::graphics::plot(Vector<double> &x, Vector<double> &y, Vector<double> &z, unsigned char options, int argc, const char *argv[])
{
    OpenGL::setData(x, y, z);
    OpenGL_Plot._options = options;
    
    glutInit(&argc, (char **)argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);
    glutInitWindowSize(OpenGL_Plot._W, OpenGL_Plot._H);
    glutInitWindowPosition(100, 100);
    glutCreateWindow((char *)OpenGL_Plot._windowName.c_str());
    OpenGL::initRendering();
    glutDisplayFunc(OpenGL::display3D);
    glutTimerFunc(OpenGL_Plot._dt, OpenGL::timerFunction2D, 0);
    glutTimerFunc(OpenGL_Plot._dt, OpenGL::refreshCameraPosition, 0);
    glutKeyboardFunc(OpenGL::handleKeyPress);
    glutKeyboardUpFunc(OpenGL::handleKeyRelease);
    glutSpecialFunc(OpenGL::handleSpecialKeyPress);
    glutSpecialUpFunc(OpenGL::handleSpecialKeyReleased);
    glutReshapeFunc(OpenGL::reshape3D);
    glutMainLoop();
}

void cda::graphics::plot(const Vector<double> &x, const Vector<double> &y, const Matrix<double> &z, unsigned char options, int argc, const char *argv[])
{
    OpenGL::setData(x, y, z);
    OpenGL_Plot._options = options;
    
    glutInit(&argc, (char **)argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);
    glutInitWindowSize(OpenGL_Plot._W, OpenGL_Plot._H);
    glutInitWindowPosition(100, 100);
    glutCreateWindow((char *)OpenGL_Plot._windowName.c_str());
    OpenGL::initRendering();
    glutDisplayFunc(OpenGL::display3D);
    glutTimerFunc(OpenGL_Plot._dt, OpenGL::timerFunction3D, 1);
    glutTimerFunc(OpenGL_Plot._dt, OpenGL::refreshCameraPosition, 2);
    glutKeyboardFunc(OpenGL::handleKeyPress);
    glutKeyboardUpFunc(OpenGL::handleKeyRelease);
    glutSpecialFunc(OpenGL::handleSpecialKeyPress);
    glutSpecialUpFunc(OpenGL::handleSpecialKeyReleased);
    glutReshapeFunc(OpenGL::reshape3D);
    glutMainLoop();
}
