//
//  MyOpenGL.h
//  Vibration Normal Modes 2D
//
//  Created by Carlos David on 5/8/12.
//  Copyright (c) 2012 NEPTUNO. All rights reserved.
//

#pragma once

#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>

#import "Vectorial.cpp"

#define OpenGL_Class OpenGL::MyOpenGL_Class
#define ColorMap_Class OpenGL::ColorMap

//  OPENGL
namespace OpenGL {
    
    //  OPENGL - CLASS
    
#define Points      0x01
#define Lines       0x02
#define Shading     0x04
#define Axis_On     0x08

#define Rot_D       0x01
#define Zoom_D      0x02
#define WN_D        0x04
#define dt_D        0x08
#define WSize_D     0x10
#define Update_D    0x20
#define SetDefault  0x40
    
    class MyOpenGL_Class {
    private:
        float _rotX, _rotY, _rotZ;
        float _zoomX, _zoomY, _zoomZ;
        float _up_down, _left_right;
        float _limXmin, _limYmin, _limZmin, _limXmax, _limYmax, _limZmax;
        float _cX, _cY, _cZ;
        float _red, _green, _blue;
        float _vMax, _vMin;
        vector<double> _vX, _vY, _vZ;
        matrix<double> _mZ;
        
        void (* _updateData)();
        
    public:
        //  Plot options
        unsigned char _options;
        string _windowName;
        int _dt, _plots_calc;
        float _W, _H;
        
        //  Init class
        static void initRendering();
        
        
        //  Set params
        static void setData(vector<double> &x, vector<double> &y, unsigned char setDefault);
        static void setData(vector<double> &x, vector<double> &y);
        static void setData(vector<double> &x, vector<double> &y, vector<double> &z);
        static void setData(vector<double> &x, vector<double> &y, vector<double> &z, unsigned char setDefault);
        static void setData(vector<double> &x, vector<double> &y, matrix<double> &z);
        static void setData(vector<double> &x, vector<double> &y, matrix<double> &z, unsigned char setDefault);
        
        //  Set adjustments
        static void setColors(float red, float green, float blue);
        
        static void setZoomsD(float zoomX, float zoomY, unsigned char setDefault);
        static void setZooms(float zoomX, float zoomY);
        static void setZoomsD(float zoomX, float zoomY, float zoomZ, unsigned char setDefault);
        static void setZooms(float zoomX, float zoomY, float zoomZ);
        static void setAxisDirection(float rotX, float rotY, float rotZ, unsigned char setDefault);
        static void setAxisDirection(float rotX, float rotY, float rotZ);
        static void setStepTime(float dt_milli, unsigned char setDefault);
        static void setStepTime(float dt_milli);
        static void setWindowSize(float Width, float Height, unsigned char setDefault);
        static void setWindowSize(float Width, float Height);
        static void setUpdateData(void updateData(), unsigned char setDefault);
        static void setUpdateData(void updateData());
        static void setWindowName(string windowName, unsigned char setDefault);
        static void setWindowName(string windowName);
        
        static void setDefaultParameters(unsigned char params);
        
        
        //  Plot functions
        static void reshape2D(int Width, int Height);
        static void reshape3D(int Width, int Height);
        static void refreshCameraPosition(int value);
        
        //  -- 2D functions --
        static void display2D();
        static void drawSolution2D();
        static void timerFunction2D(int value);
        static void plotAxes2D();
        static void center2D();
        
        //  -- 3D functions
        static void display3D();
        static void drawSolution3D();
        static void timerFunction3D(int value);
        static void plotAxes3D();
        static void center3D();
        
        
        //  Keyboard and mouse functions
#define UP              0x01
#define DOWN            0x02
#define RIGHT           0x04
#define LEFT            0x08
        
#define Z_in            0x10
#define Z_out           0x20
#define A_              0x40
#define D_              0x80
        
        static void handleKeyPress(unsigned char key, int x, int y);
        static void handleKeyRelease(unsigned char key, int x, int y);
        static void handleSpecialKeyPress(int key, int x, int y);
        static void handleSpecialKeyReleased(int key, int x, int y);
        
    };
    
    //  COLORMAPS FOR OPENGL - CLASS
    
#define JET             0x01
#define HOT_to_COLD     0x02
#define COOL            0x04
#define WINTER          0x08
#define HOT             0x10
#define SUMMER          0x20
#define JET_FIXED       0x40
    
    class ColorMap {
    private:
        static float red, green, blue;
        
        static void Jet(float value, float min, float max);
        static void Hot_to_Cold(float value, float min, float max);
        static void Cool(float value, float min, float max);
        static void Winter(float value, float min, float max);
        static void Hot(float value, float min, float max);
        static void Summer(float value, float min, float max);
        static void Jet_Fixed(float value, float min, float max);
    public:
        static unsigned char colormap;
        static void setColor(float value, float min, float max);
    };
    
    //  PLOT FUNCTIONS - NAMESPACE
    void plot(vector<double> &x, vector<double> &y, unsigned char options, int argc, const char *argv[]);
    void plot(vector<double> &x, vector<double> &y, vector<double> &z, unsigned char options, int argc, const char *argv[]);
    void plot(vector<double> &x, vector<double> &y, matrix<double> &z, unsigned char options, int argc, const char *argv[]);
}
