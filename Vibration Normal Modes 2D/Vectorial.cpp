//
//  Vectorial.cpp
//  Class for vectorial operations
//  Vibration Normal Modes 2D
//
//  Created by Carlos David on 5/8/12.
//  Copyright (c) 2012 NEPTUNO. All rights reserved.
//

#include "Vectorial.h"

#include <typeinfo>


const int tolMax = 1000;
const double preDef = 1E-04;
const string RV = "\n[vector::";
const string RM = "\n[matrix::";

                //  --- VECTOR CLASS ---

//  --- DEFINITION ---
MAT_TEMPLATE inline
vectorT::vector(int _n)
{
    n = _n;
    v = new T[n];
}

MAT_TEMPLATE inline
vectorT::vector()
{
    n = 1;
    v = new T[n];
}

MAT_TEMPLATE inline
vectorT::vector(T* _v)
{
    v = _v;
}

MAT_TEMPLATE inline
vectorT::vector(const vectorT& V)
{
    n = V.n;
    v = new T[n];
    copy(V.v, V.v + V.n, v);
}

MAT_TEMPLATE inline
vectorT::~vector()
{
    delete [] v;
}

MAT_TEMPLATE inline void
vectorT::reSize(int _n)
{
    if (n != _n) {
        vectorT tmp(n);
        tmp = (* this);
        
        delete [] v;
        n = _n;
        v = new T[n];
        
        this -> set(0, tmp, n);
        
        for (int i=tmp.n; i<_n; i++) {
            v[i] = (T)0.0;
        }
    }
}

MAT_TEMPLATE inline vectorT&
vectorT::operator=(const vectorT& V)
{
    if (&V != this) {
        delete[] v;
        n = V.n;
        v = new T[n];
        copy(V.v, V.v + V.n, v);
    }
    return (* this);
}



//  --- GETS AND SETS ---
//  Gets
MAT_TEMPLATE inline vectorT
vectorT::get(int _i, int lenght)
{
    if (n-_i < lenght)
        cout << RV << "getElements(i, lenght)] - El vector que se intenta copiar no tiene tantos elementos.\n";
    
    vectorT tmp(n-_i);
    tmp.zero();
    
    if (n-_i > lenght)
        for (int i=_i; i<lenght; i++) {
            tmp.v[i-_i] = v[i];
        }
    else
        for (int i=_i; i<n; i++) {
            tmp.v[i-_i] = v[i];
        }
    
    return tmp;
}

MAT_TEMPLATE inline vectorT
vectorT::get(int _i)
{
    return (this -> get(_i, n-_i));
}

//  Sets
MAT_TEMPLATE inline void
vectorT::set(int _i, const vectorT& V, int lenght)
{
    if (lenght > n) {
        for (int i=_i; i<n; i++) {
            v[i] = V.v[i-_i];
        }
        cout << RV << "setElements(i, V)] - El vector que se ha intentado añadir es demasiado grande. No se ha añadido todo el vector.\n";
    } else {
        for (int i=_i; i<lenght+_i; i++) {
            v[i] = V.v[i-_i];
        }
    }
}

MAT_TEMPLATE inline void
vectorT::set(int _i, const vectorT& V)
{
    this -> set(_i, V, V.n);
}



//  --- FUNCTIONS ---
//  Returns an integer / a float / a double
MAT_TEMPLATE inline int
vectorT::dim() const
{
    return n;
}

MAT_TEMPLATE inline T
vectorT::max() const
{
	T max = v[0];
	for (int i=1; i<n; i++) {
		if (v[i] > max) {
			max = v[i];
		}
	}
	return max;
}

MAT_TEMPLATE inline T
vectorT::min() const
{
    T min = v[0];
	for (int i=1; i<n; i++) {
		if (v[i] < min) {
			min = v[i];
		}
	}
	return min;
}

MAT_TEMPLATE inline T
vectorT::maxAbs() const
{
    T maxA = fabs(v[0]);
    for (int i=1; i<n; i++) {
        if (fabs(v[i]) > maxA) {
            maxA = fabs(v[i]);
        }
    }
    return maxA;
}

MAT_TEMPLATE inline T
vectorT::sum() const
{
    T sum = (T)0.0;
    for (int i=0; i<n; i++) {
        sum += v[i];
    }
    return sum;
}

MAT_TEMPLATE inline T&
vectorT::operator()(int _i) const
{
    return v[_i];
}


//  Returns a vector / a matrix
MAT_TEMPLATE inline vectorT
vectorT::mod(const T& D)
{
    T norm = !(* this);
    for (int i=0; i<n; i++) {
        v[i] *= D / norm;
    }
    return (* this);
}

MAT_TEMPLATE inline vectorT
vectorT::unitary()
{
    T sum = (T)0.0;
    for (int i=0; i<n; i++) {
        sum += v[i]*v[i];
    }
    return ((* this)/std::sqrt(sum));
}

MAT_TEMPLATE inline vectorT
vectorT::organize()
{
    vectorT tmp(n);
    bool change = true;
    T vb = (T)0.0;
    tmp = (* this);
    
    while (change) {
        change = false;
        for (int i=1; i<n; i++) {
            if (tmp.v[i-1] > tmp.v[i]) {
                vb = tmp.v[i-1];
                tmp.v[i-1] = tmp.v[i];
                tmp.v[i] = vb;
                change = true;
            }
        }
    }
    
    return tmp;
}

MAT_TEMPLATE inline matrixT
vectorT::transpose()
{
    matrixT tmp(n,1);
    tmp.zero();
    for (int i=0; i<n; i++) {
        tmp(i,0) = v[i];
    }
    return tmp;
}

MAT_TEMPLATE inline vectorT
vectorT::zero(int _n)
{
    vectorT tmp(_n);
    tmp.zero();
    return tmp;
}

MAT_TEMPLATE inline vectorT
vectorT::ones(int _n)
{
    vectorT tmp(_n);
    tmp.ones();
    return tmp;
}


//  Returns a bool
MAT_TEMPLATE inline bool
vectorT::null()
{
    for (int i=0; i<n; i++) {
        if (v[i] != (T)0.0)
            return false;
    }
    return true;
}

MAT_TEMPLATE inline bool
vectorT::duplicate(T range)
{
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            if (i != j && (v[j] >= v[i]*(1.0-range) && v[j] <= v[i]*(1.0+range)))
                return true;
        }
    }
    return false;
}

MAT_TEMPLATE inline bool
vectorT::duplicate()
{
    return (this -> duplicate(0.0));
}


//  Void functions
MAT_TEMPLATE inline void
vectorT::ones()
{
    for (int i=0; i<n; i++) {
        v[i] = (T)1.0;
    }
}

MAT_TEMPLATE inline void
vectorT::zero()
{
    for (int i=0; i<n; i++) {
        v[i] = (T)0.0;
    }
}

MAT_TEMPLATE inline void
vectorT::rand()
{
    for (int i=0; i<n; i++) {
        v[i] = std::rand();
    }
}

MAT_TEMPLATE inline void
vectorT::rand(const T& Min, const T& Max)
{
    for (int i=0; i<n; i++) {
        v[i] = drand48() * (Max - Min) + Min;
    }
}

MAT_TEMPLATE inline void
vectorT::write(const string& path, const string& filename)
{
    const string Npath = path + filename;
    ofstream out(Npath.data());
    out << (* this);
}

MAT_TEMPLATE inline void
vectorT::write(const string& filename)
{
    const string home = getenv("HOME");
    const string path = home + "/Desktop/";
    write(path, filename);
}



//  --- OPERATORS ---
MAT_TEMPLATE inline vectorT
vectorT::operator+(const vectorT& V)
{
    vectorT tmp(n);
    tmp.zero();
    if (n != V.n) {
        cout << RV << "operator+V] - Las dimensiones de los vectores no coinciden.\n";
        tmp.ones();
    } else {
        for (int i=0; i<n; i++) {
            tmp.v[i] = v[i] + V.v[i];
        }
    }
    return tmp;
}

MAT_TEMPLATE inline vectorT&
vectorT::operator+=(const vectorT& V)
{
    if (n != V.n) {
        cout << RV << "operator+=V] - Las dimensiones de los vectores no coinciden.\n";
    }
    for (int i=0; i<n; i++) {
        v[i] += V.v[i];
    }
    return (* this);
}

MAT_TEMPLATE inline vectorT
vectorT::operator-(const vectorT& V)
{
    vectorT tmp(n);
    tmp.zero();
    if (n != V.n) {
        cout << RV << "operator-V] - Las dimensiones de los vectores no coinciden.\n";
        tmp.ones();
    } else {
        for (int i=0; i<n; i++) {
            tmp.v[i] = v[i] - V.v[i];
        }
    }
    return tmp;
}

MAT_TEMPLATE inline vectorT&
vectorT::operator-=(const vectorT& V)
{
    if (n != V.n) {
        cout << RV << "operator-=V] - Las dimensiones de los vectores no coinciden.\n";
    }
    for (int i=0; i<n; i++) {
        v[i] -= V.v[i];
    }
    return (* this);
}

MAT_TEMPLATE inline vectorT
vectorT::operator*(const T& D)
{
    vectorT tmp(n);
    tmp.zero();
    for (int i=0; i<n; i++) {
        tmp.v[i] = v[i] * D;
    }
    return tmp;
}

MAT_TEMPLATE inline vectorT&
vectorT::operator*=(const T& D)
{
    for (int i=0; i<n; i++) {
        v[i] *= D;
    }
    return (* this);
}

MAT_TEMPLATE inline vectorT
vectorT::operator*(const matrixT& Mc)
{
    int m = Mc.columns();
    vectorT tmp(m);
    tmp.zero();
    if (n != Mc.rows()) {
        cout << RV << "operator*Mc] - La altura de la matriz no corresponde con la longitud del vector.\n";
    } else {
        for (int j=0; j<m; j++) {
            for (int k=0; k<n; k++) {
                tmp.v[j] += v[k]*Mc(k,j);
            }
        }
    }
    return tmp;
}

MAT_TEMPLATE inline vectorT
vectorT::operator/(const T& D)
{
    vectorT tmp(n);
    tmp.zero();
    for (int i=0; i<n; i++) {
        tmp.v[i] = v[i]/D;
    }
    return tmp;
}

MAT_TEMPLATE inline vectorT&
vectorT::operator/=(const T& D)
{
    for (int i=0; i<n; i++) {
        v[i] /= D;
    }
    return (* this);
}

MAT_TEMPLATE inline vectorT
vectorT::operator%(const int& I)
{
    vectorT tmp(n);
    for (int i=0; i<n; i++) {
        tmp.v[i] = ((int)v[i])%I;
    }
    return tmp;
}

MAT_TEMPLATE inline vectorT&
vectorT::operator%=(const int& I)
{
    for (int i=0; i<n; i++) {
        v[i] = ((int)v[i])%I;
    }
    return (* this);
}

MAT_TEMPLATE inline vectorT
vectorT::operator-()
{
    vectorT tmp(n);
    tmp.zero();
    for (int i=0; i<n; i++) {
        tmp.v[i] = -v[i];
    }
    return tmp;
}

MAT_TEMPLATE inline vectorT
vectorT::cross(const vectorT& V)
{
    vectorT tmp(n);
    tmp.zero();
    if (n == 3 && V.n == 3) {
        tmp.v[0] = v[1]*V.v[2] - v[2]*V.v[1];
        tmp.v[1] = v[2]*V.v[0] - v[0]*V.v[2];
        tmp.v[2] = v[0]*V.v[1] - v[1]*V.v[0];
    } else {
        cout << RV << "cross(V)] - Los vectores no son tridimensionales.\n";
        tmp.ones();
    }
    return tmp;
}

MAT_TEMPLATE inline T
vectorT::operator!()
{
    return std::sqrt(~(* this));
}

MAT_TEMPLATE inline T
vectorT::operator~()
{
    return ((* this) * (* this));
}


//  --- MORE OPERATORS ---
MAT_TEMPLATE inline T
vectorT::operator*(const vectorT& V)
{
    T tmp = (T)0.0;
    if (n != V.n) {
        cout << RV << "operator*V] - Las dimensiones de los vectores no coinciden.\n";
    } else {
        for (int i=0; i<n; i++) {
            tmp += v[i] * V.v[i];
        }
    }
    return tmp;
}

MAT_TEMPLATE inline vectorT
operator*(const T& D, const vectorT& V)
{
    vectorT tmp(V.dim());
    tmp.zero();
    for (int i=0; i<V.dim(); i++) {
        tmp(i) = V(i) * D;
    }
    return tmp;
}

MAT_TEMPLATE inline void
operator>>(istream& in, vectorT& V)
{
    int n=1, k=0, cDim;
    string c, num;
    
    getline(in, c);
    cDim = (unsigned int)c.size();
    for (int i=0; i<cDim; i++) {
        if (c[i] == ',')
            n++;
    }
    
    V.reSize(n);
    for (int j=0; j<n; j++) {
        num = "";
        while (c[k] != ',' && (c[k] != '\n' && c[k] != '\r')) {
            num += c[k];
            k++;
            if (k == cDim)
                break;
        }
        k++;
        V(j) = strtod(num.c_str(), NULL);
    }
}

MAT_TEMPLATE inline ostream&
operator<<(ostream& out, const vectorT& V)
{
//    if (out == cout) {
//        out.precision(6);
//        out << "[" << V(0);
//        if (V.dim() > 1) {
//            for (int i=1; i<V.dim(); i++) {
//                out << ", " << V(i); 
//            }
//        }
//        out << "]\n";
//    } else {
//        out << V(0);
//        if (V.dim() > 1) {
//            for (int i=1; i<V.dim(); i++) {
//                out << "," << V(i); 
//            }
//        }
//    }
    return out;
}



//  --- OTHER FUNCTIONS ---
MAT_TEMPLATE inline vectorT
sqrt(const vectorT& V)
{
    vectorT tmp(V.dim());
    for (int i=0; i<V.dim(); i++) {
        tmp(i) = std::sqrt(V(i));
    }
    return tmp;
}

MAT_TEMPLATE inline vectorT
pow(const vectorT& V, double exp)
{
    vectorT tmp(V.dim());
    for (int i=0; i<V.dim(); i++) {
        tmp(i) = std::pow(V(i), exp);
    }
    return tmp;
}

MAT_TEMPLATE inline vectorT
log(const vectorT& V)
{
    vectorT tmp(V.dim());
    for (int i=0; i<V.dim(); i++) {
        tmp(i) = std::log(V(i));
    }
    return tmp;
}

MAT_TEMPLATE inline vectorT
zero(int _n)
{
    vectorT tmp(_n);
    tmp.zero();
    return tmp;
}

MAT_TEMPLATE inline vectorT
ones(int _n)
{
    vectorT tmp(_n);
    tmp.ones();
    return tmp;
}

MAT_TEMPLATE inline vectorT
rand(int _n)
{
    vectorT tmp(_n);
    tmp.rand();
    return tmp;
}

MAT_TEMPLATE inline vectorT
rand(int _n, const T& Min, const T& Max)
{
    vectorT tmp(_n);
    tmp.rand(Min, Max);
    return tmp;
}

MAT_TEMPLATE inline vectorT
round(const vectorT& V)
{
    vectorT tmp(V.dim());
    tmp.zero();
    for (int i=0; i<V.dim(); i++) {
        if (V(i) >= 0 && V(i)-floor(V(i)) < 0.5)
            tmp(i) = (int)floor(V(i));
        else if (V(i) >= 0 && V(i)-floor(V(i)) >= 0.5)
            tmp(i) = (int)ceil(V(i));
        else if (V(i) < 0 && V(i)-floor(V(i)) >= 0.5)
            tmp(i) = (int)ceil(V(i));
        else
            tmp(i) = (int)floor(V(i));
    }
    
    return tmp;
}



            //  --- MATRIX CLASS ---

//  --- DEFINITION ---
MAT_TEMPLATE inline
matrixT::matrix(int _n, int _m)
{
    n = _n;
    m = _m;
    size = n*m;
    a = new T[size];
}

MAT_TEMPLATE inline
matrixT::matrix()
{
    n = 1;
    m = 1;
    size = n*m;
    a = new T[size];
}

MAT_TEMPLATE inline
matrixT::matrix(T* _a)
{
    a = _a;
}

MAT_TEMPLATE inline
matrixT::matrix(const matrixT& M)
{
    n = M.n;
    m = M.m;
    size = M.size;
    a = new T[size];
    copy(M.a, M.a + M.size, a);
}

MAT_TEMPLATE
template <class T2>
matrixT::matrix(const matrix<T2>& M)
{
    n = M.rows();
    m = M.columns();
    size = M.elements();
    a = new T[size];
    
    for (int k=0; k<size; k++) {
        a[k] = (T)M(k);
    }
}

MAT_TEMPLATE inline
matrixT::~matrix()
{
    delete [] a;
}

MAT_TEMPLATE inline void
matrixT::reSize(int _n, int _m)
{
    matrixT tmp(n,m);
    tmp = (* this);
    
    delete [] a;
    n = _n;
    m = _m;
    size = n*m;
    a = new T[size];
    
    this -> setMatrix(0, 0, tmp);
    
    for (int i=0; i<tmp.n; i++) {
        for (int j=tmp.m; j<_m; j++) {
            a[i*_m + j] = (T)0.0;
        }
    }
    
    for (int i=tmp.n; i<_n; i++) {
        for (int j=0; j<_m; j++) {
            a[i*_m + j] = (T)0.0;
        }
    }
}

MAT_TEMPLATE inline matrixT&
matrixT::operator=(const matrixT& M)
{
    if (&M != this) {
        delete [] a;
        n = M.n;
        m = M.m;
        size = M.size;
        a = new T[size];
        copy(M.a, M.a + M.size, a);
    }
    return (* this);
}



//  --- GETS AND SETS ---
//  Gets
MAT_TEMPLATE inline matrixT
matrixT::getRow(int _i)
{
    matrixT tmp(1,m);
    tmp.zero();
    if (_i >= n) {
        cout << RM << "getRow(i)] - Esta fila no existe, se ha devuelvo una fila de unos.\n";
        tmp.ones();
    } else {
        for (int j=0; j<m; j++) {
            tmp.a[j] = a[_i*m + j];
        }
    }
    return tmp;
}

MAT_TEMPLATE inline matrixT
matrixT::getColumn(int _j)
{
    matrixT tmp(n,1);
    tmp.zero();
    if (_j >= m) {
        cout << RM << "getColumn(j)] - Esta columna no existe, se ha devuelvo una columna de unos.\n";
        tmp.ones();
    } else {
        for (int i=0; i<n; i++) {
            tmp.a[i] = a[i*m + _j];
        }
    }
    return tmp;
}

MAT_TEMPLATE inline matrixT
matrixT::getMatrix(int _i, int _j, int height, int columns)
{
    matrixT tmp(height,columns);
    tmp.zero();
    if (n < _i+height || m < _j+columns) {
        cout << RM << "getMatrix(i, j, height, columns)] - Estás intentado acceder a una parte de la matriz que no existe, se ha devuelto una matriz de unos\n";
        tmp.ones();
    } else {
        for (int i=0; i<height; i++) {
            for (int j=0; j<columns; j++) {
                tmp.a[i*columns + j] = a[(i+_i)*m + (j+_j)];
            }
        }
    }
    return tmp;
}

MAT_TEMPLATE inline matrixT
matrixT::getMatrix(int _i, int _j)
{
    return (this -> getMatrix(_i, _j, n-_i, m-_j));
}

MAT_TEMPLATE inline vectorT
matrixT::getDiagonal()
{
    vectorT tmp(n);
    tmp.zero();
    if (n != m)
        cout << RM << "getDiagonal()] - Esta no es una matriz cuadrada.\n";
    for (int i=0; i<n; i++) {
        tmp(i) = a[i*m + i];
    }
    return tmp;
}

MAT_TEMPLATE inline vectorT
matrixT::getRowV(int _i, int _j)
{
    vectorT tmp(m-_j);
    tmp.zero();
    if (_i >= n) {
        cout << RM << "getRowV(i, j)] - Esta fila no existe, se ha devuelvo una fila de unos.\n";
        tmp.ones();
    } else {
        for (int j=0; j<m-_j; j++) {
            tmp(j) = a[_i*m + (j+_j)];
        }
    }
    return tmp;
}

MAT_TEMPLATE inline vectorT
matrixT::getRowV(int _i)
{
    return (this -> getRowV(_i, 0));
}

MAT_TEMPLATE inline vectorT
matrixT::getColumnV(int _i, int _j)
{
    vectorT tmp(n-_i);
    tmp.zero();
    if (_j >= m) {
        cout << RM << "getColumnV(i, j)] - Esta columna no existe, se ha devuelvo una columna de unos.\n";
        tmp.ones();
    } else {
        for (int i=0; i<n-_i; i++) {
            tmp(i) = a[(i+_i)*m + _j];
        }
    }
    return tmp;
}

MAT_TEMPLATE inline vectorT
matrixT::getColumnV(int _j)
{
    return (this -> getColumnV(0, _j));
}


//  Sets
MAT_TEMPLATE inline void
matrixT::setColumn(int _i, int _j, const matrixT& Mc)
{
    if (n < _i+Mc.n || Mc.m > 1)
        cout << RM << "setColumn(i, j, Mc)] - La columna que se está intentando sustituir es demasiado larga o no se está introduciendo una matriz columna.\n";
    for (int i=_i; i<_i+Mc.n; i++) {
        a[i*m + _j] = Mc.a[(i-_i)*Mc.m];
    }
}

MAT_TEMPLATE inline void
matrixT::setColumn(int _j, const matrixT& Mc)
{
    this -> setColumn(0, _j, Mc);
}

MAT_TEMPLATE inline void
matrixT::setColumnV(int _i, int _j, const vectorT& V)
{
    if (n < _i+V.dim())
        cout << RM << "setColumnV(i, j, V)] - La columna que se está intentando sustituir es demasiado larga.\n";
    for (int i=_i; i<_i+V.dim(); i++) {
        a[i*m + _j] = V(i-_i);
    }
}

MAT_TEMPLATE inline void
matrixT::setColumnV(int _j, const vectorT& V)
{
    this -> setColumnV(0, _j, V);
}

MAT_TEMPLATE inline void
matrixT::setRow(int _i, int _j, const matrixT& Mr)
{
    if (m < _j+Mr.m || Mr.n > 1) {
        cout << RM << "setRow(i, j, Mr)] - La fila que se está intentando sustituir es demasiado larga o no se está introduciendo una matriz fila.\n";
    }
    for (int j=_j; j<_j+Mr.m; j++) {
        a[_i*m + j] = Mr.a[j-_j];
    }
}

MAT_TEMPLATE inline void
matrixT::setRow(int _i, const matrixT& Mr)
{
    this -> setRow(_i, 0, Mr);
}

MAT_TEMPLATE inline void
matrixT::setRowV(int _i, int _j, const vectorT& V)
{
    if (m < _j+V.dim())
        cout << RM << "setRowV(i, j, V)] - La fila que se está intentando sustituir es demasiado larga.\n";
    for (int j=_j; j<_j+V.dim(); j++) {
        a[_i*m + j] = V(j-_j);
    }
}

MAT_TEMPLATE inline void
matrixT::setRowV(int _i, const vectorT& V)
{
    this -> setRowV(_i, 0, V);
}

MAT_TEMPLATE inline void
matrixT::setMatrix(int _i, int _j, const matrixT& M)
{
    if (n < M.n+_i && m < M.m+_j)
        for (int i=_i; i<n; i++) {
            for (int j=_j; j<m; j++) {
                a[i*m + j] = M.a[(i-_i)*M.m + (j-_j)];
            }
        }
    
    else if (n < M.n)
        for (int i=_i; i<n; i++) {
            for (int j=_j; j<M.m+_j; j++) {
                a[i*m + j] = M.a[(i-_i)*M.m + (j-_j)];
            }
        }
    
    else if (m < M.m+_j)
        for (int i=_i; i<M.n+_i; i++) {
            for (int j=_j; j<m; j++) {
                a[i*m + j] = M.a[(i-_i)*M.m + (j-_j)];
            }
        }
    
    else
        for (int i=_i; i<M.n+_i; i++) {
            for (int j=_j; j<M.m+_j; j++) {
                a[i*m + j] = M.a[(i-_i)*M.m + (j-_j)];
            }
        }
}



//  --- FUNCTIONS ---
//  Returns an integer / a double / a float
MAT_TEMPLATE inline int
matrixT::dims() const
{
    int dims[2];
    dims[0] = n;
    dims[1] = m;
    return (* dims);
}

MAT_TEMPLATE inline int
matrixT::rows() const
{
    return n;
}

MAT_TEMPLATE inline int
matrixT::columns() const
{
    return m;
}

MAT_TEMPLATE inline int
matrixT::elements() const
{
    return size;
}


MAT_TEMPLATE inline T
matrixT::sumRow(int _i) const
{
    if (n >= _i)
        cout << RM << "sumRow(i)] - Esta fila no existe.\n";
    T sum = (T)0.0;
    for (int j=0; j<m; j++) {
        sum += a[_i*m + j];
    }
    return sum;
}

MAT_TEMPLATE inline T
matrixT::sumColumn(int _j) const
{
    if (m >= _j)
        cout << RM << "sumColumn(j)] - Esta columna no existe.\n";
    T sum = (T)0.0;
    for (int i=0; i<n; i++) {
        sum += a[i*m + _j];
    }
    return sum;
}

MAT_TEMPLATE inline T
matrixT::max() const
{
    T max = a[0];
    for (int k=1; k<size; k++) {
        if (max < a[k])
            max = a[k];
    }
    return max;
}

MAT_TEMPLATE inline T
matrixT::maxAbs() const
{
    T max = fabs(a[0]);
    for (int k=1; k<size; k++) {
        if (fabs(a[k]) > max)
            max = fabs(a[k]);
    }
    return max;
}

MAT_TEMPLATE inline T
matrixT::maxAbs_sig() const
{
    T max = a[0];
    for (int k=1; k<size; k++) {
        if (fabs(a[k]) > fabs(max))
            max = a[k];
    }
    return max;
}

MAT_TEMPLATE inline T
matrixT::min() const
{
    T min = a[0];
    for (int k=1; k<size; k++) {
        if (min > a[k])
            min = a[k];
    }
    return min;
}

MAT_TEMPLATE inline T
matrixT::det()
{
    T det = 1.0;
    if (n != m) {
        cout << RM << "det()] - La matriz no es cuadrada, así que no se puede calcular su determinante.\n";
        return det;
    }
    matrixT tmp(n,m);
    tmp.zero();
    tmp = this -> U();
    
    for (int i=0; i<n; i++) {
        det *= tmp(i,i);
    }
    
    return det;
}

MAT_TEMPLATE inline T&
matrixT::operator()(int _i, int _j) const
{
    if (n <= _i || m <= _j)
        cout << RM << "operator(i,j)] - Este elemento de matriz no existe.\n";
    return a[_i*m + _j];
}

MAT_TEMPLATE inline T&
matrixT::operator()(int _k) const
{
    if (size <= _k)
        cout << RM << "operator(k)] - Este elemento de matriz no existe.\n";
    return a[_k];
}


//  Returns a matrix
MAT_TEMPLATE inline matrixT
matrixT::sumRows()
{
    T sum;
    matrixT tmp(n,1);
    tmp.zero();
    for (int i=0; i<n; i++) {
        sum = (T)0.0;
        for (int j=0; j<m; j++) {
            sum += a[i*m + j];
        }
        tmp.a[i] = sum;
    }
    return tmp;
}

MAT_TEMPLATE inline matrixT
matrixT::sumColumns()
{
    T sum;
    matrixT tmp(1,m);
    tmp.zero();
    for (int i=0; i<n; i++) {
        sum = (T)0.0;
        for (int j=0; j<m; j++) {
            sum += a[j*m + i];
        }
        tmp.a[i] = sum;
    }
    return tmp;
}

MAT_TEMPLATE inline matrixT
matrixT::transpose()
{
    matrixT tmp(m,n);
    tmp.zero();
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            tmp.a[i*n + j] = a[j*m + i];
        }
    }
    return tmp;
}

MAT_TEMPLATE inline matrixT
matrixT::zero(int _n, int _m)
{
    matrixT tmp(_n,_m);
    tmp.zero();
    return tmp;
}

MAT_TEMPLATE inline matrixT
matrixT::ones(int _n, int _m)
{
    matrixT tmp(_n,_m);
    tmp.ones();
    return tmp;
}

MAT_TEMPLATE inline matrixT
matrixT::identity(int _n, int _m)
{
    if (_n != _m)
        cout << RM << "identity(n, m)] - No es una matriz cuadrada.\n";
    matrixT I(_n,_m);
    I.identity();
    return I;
}


MAT_TEMPLATE inline vectorT
matrixT::sumRowsV()
{
    T sum;
    vectorT tmp(n);
    tmp.zero();
    for (int i=0; i<n; i++) {
        sum = (T)0.0;
        for (int j=0; j<m; j++) {
            sum += a[i*m + j];
        }
        tmp(i) = sum;
    }
    return tmp;
}

MAT_TEMPLATE inline vectorT
matrixT::sumColumnsV()
{
    T sum;
    vectorT tmp(m);
    tmp.zero();
    for (int i=0; i<n; i++) {
        sum = (T)0.0;
        for (int j=0; j<m; j++) {
            sum += a[j*m + i];
        }
        tmp(i) = sum;
    }
    return tmp;
}


//  Returns a bool
MAT_TEMPLATE inline bool
matrixT::null()
{
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            if (a[i*m + j] != (T)0.0)
                return false;
        }
    }
    return true;
}

MAT_TEMPLATE inline bool
matrixT::duplicate(T range)
{
    for (int k=0; k<size; k++) {
        for (int p=0; p<size; p++) {
            if (k != p && (a[p] >= a[k]*(1.0-range) && a[p] <= a[k]*(1.0+range)))
                return true;
        }
    }
    return false;
}

MAT_TEMPLATE inline bool
matrixT::duplicate()
{
    return (this -> duplicate(0.0));
}


//  Void functions
MAT_TEMPLATE inline void
matrixT::zero()
{
    for (int k=0; k<size; k++) {
        if (typeid(T) == typeid(bool))
            a[k] = false;
        else
            a[k] = (T)0.0;
    }
}

MAT_TEMPLATE inline void
matrixT::ones()
{
    for (int k=0; k<size; k++) {
        if (typeid(T) == typeid(bool))
            a[k] = true;
        else
            a[k] = (T)1.0;
    }
}

MAT_TEMPLATE inline void
matrixT::identity()
{
    if (n != m)
        cout << RM << "identity()] - No es una matriz cuadrada.\n";
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            if (i == j) {
                a[i*m + j] = (T)1.0;
            } else {
                a[i*m + j] = (T)0.0;
            }
        }
    }
}

MAT_TEMPLATE inline void
matrixT::write(const string& path, const string& filename)
{
    const string Npath = path + filename;
    ofstream out(Npath.data());
    out << (* this);
}

MAT_TEMPLATE inline void
matrixT::write(const string& filename)
{
    const string home = getenv("HOME");
    const string path = home + "/Desktop/";
    write(path, filename);
}



//  --- MÉTODO LU ---
MAT_TEMPLATE inline matrixT
matrixT::LU()
{
    if (n != m) {
        cout << RM << "LU()] - No se puede calcular la diagonal inferior porque no es una matriz cuadrada.\n";
        matrixT tmp(n,m);
        tmp.ones();
        return tmp;
    }
    
    matrixT LU(n,m);
    LU.zero();
    LU = *this;
    
    for (int i=1; i<n; i++) {
        LU.a[i*m] /= LU.a[0];
    }
    T sum = (T)0.0;
    for (int i=1; i<n-1; i++) {
        for (int j=i; j<m; j++) {
            sum = (T)0.0;
            for (int k=0; k<i; k++) {
                sum += LU.a[k*m + j] * LU.a[i*m + k];
            }
            LU.a[i*m + j] -= sum;
        }
        for (int k=i+1; k<n; k++) {
            sum = (T)0.0;
            for (int j=0; j<i; j++) {
                sum += LU.a[j*m + i] * LU.a[k*m + j];
            }
            LU.a[k*m + i] = (LU.a[k*m + i] - sum) / LU.a[i*m + i];
        }
    }
    sum = (T)0.0;
    for (int k=0; k<n-1; k++) {
        sum += LU.a[k*m + (m-1)] * LU.a[(n-1)*m + k];
    }
    LU.a[(n-1)*m + (m-1)] -= sum;
    
    return LU;
}

MAT_TEMPLATE inline matrixT
matrixT::U()
{
    if (n != m)
        cout << RM << "U()] -> ";
    matrixT LU(n,m), U(n,m);
    LU.zero();
    U.zero();
    LU = this -> LU();
    
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            if (j >= i) {
                U.a[i*m + j] = LU.a[i*m + j];
            } else {
                U.a[i*m + j] = (T)0.0;
            }
        }
    }
    return U;
}

MAT_TEMPLATE inline matrixT
matrixT::L()
{
    if (n != m)
        cout << RM << "L()] -> ";
    matrixT LU(n,m), L(n,m);
    LU.zero();
    L.zero();
    LU = this -> LU();
    
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            if (j < i) {
                L.a[i*m + j] = LU.a[i*m + j];
            } else if (j == i) {
                L.a[i*m + j] = (T)1.0;
            } else {
                L.a[i*m + j] = (T)0.0;
            }
        }
    }
    return L;
}

MAT_TEMPLATE inline vectorT
matrixT::solveLU(const vectorT& B)
{
    vectorT tmp(n);
    tmp.zero();
    if (n != m || n != B.dim()) {
        cout << RM << "solveLU(B)] - El sistema no es compatible determinado, la matriz no es cuadrada o el vector de términos independientes no coincide con el número de incógnitas.\n";
        tmp.ones();
        return tmp;
    }
    vectorT x(n);
    x.zero();
    matrix LU(n,m);
    LU.zero();
    LU = this -> LU();
    
    T sum = (T)0.0;
    for (int i=0; i<n; i++) {
        sum = (T)0.0;
        for (int j=0; j<i; j++) {
            sum += LU.a[i*m + j] * tmp(j);
        }
        tmp(i) = B(i) - sum;
    }
    
    for (int i=n-1; i>=0; i--) {
        sum = (T)0.0;
        for (int j=m-1; j>=0; j--) {
            sum += LU.a[i*m + j] * x(j);
        }
        x(i) = (tmp(i)-sum) / LU.a[i*m + i];
    }
    return x;
}

MAT_TEMPLATE inline vectorT
matrixT::solveLU3d(const vectorT& B)
{
    vectorT tmp(n);
    tmp.zero();
    if (n != B.dim()) {
        cout << RM << "solveLU3d(B)] - El sistema no es compatible determinado, la matriz no es cuadrada o el vector de términos independientes no coincide con el número de incógnitas.\n";
        tmp.ones();
        return tmp;
    }
    vectorT al(n), be(n), x(n);
    al.zero();
    be.zero();
    x.zero();
    
    be(0) = a[1];
    for (int i=1; i<n; i++) {
        al(i) = a[i*m] / be(i-1);
        be(i) = a[i*m + 1] - al(i)*a[(i-1)*m + 2];
    }
    tmp(0) = B(0);
    for (int i=1; i<n; i++) {
        tmp(i) = B(i) - al(i)*tmp(i-1);
    }
    x(n-1) = tmp(n-1)/be(n-1);
    for (int i=(n-2); i>=0; i--) {
        x(i) = (tmp(i) - a[i*m + 2]*x(i+1)) / be(i);
    }
    
    return x;
}

MAT_TEMPLATE inline vectorT
matrixT::solveGS3d(const vectorT& B, T err)
{
    vectorT x(n);
    x.zero();
    if (n != B.dim()) {
        cout << RM << "solveGS3d(B, err)] - El sistema no es compatible determinado, la matriz no es cuadrada o el vector de términos independientes no coincide con el número de incógnitas.\n";
        x.ones();
        return x;
    }
    
    T error = 0.0, xAnt;
    do {
        error = 0.0;
        
        xAnt = x(0);
        x(0) = (B(0) - a[0]*x(1))/a[1];
        error += pow((x(0)-xAnt), 2);
        
        for (int i=1; i<(n-1); i++) {
            xAnt = x(i);
            x(i) = (B(i) - a[i*n]*x(i-1) - a[i*n + 2]*x(i+1))/a[i*n + 1];
            error += pow((x(i)-xAnt), 2);
        }
        
        xAnt = x(n-1);
        x(n-1) = (B(n-1) - a[(n-1)*n]*x(n-2))/a[(n-1)*n + 1];
        error += pow((x(n-1)-xAnt), 2);
        
        error = std::sqrt(error);
        
    } while (error > err);
    return x;
}


//  --- MÉTODO QR ---
MAT_TEMPLATE inline matrixT
matrixT::QR(unsigned char QorR)
{
    if (n != m) {
        cout << RM << "QR(QorR) - La matriz no es cuadrada.";
        return ones(n, n);
    }
    
    matrixT Q(n,n), R(n,n), I(n,n), H(n,n), A(n,n);
    Q.zero();
    R.zero();
    H.zero();
    A.zero();
    I.identity();
    
    //  First column
    vectorT c(n), vt(n);
    matrixT v(n,1);
    c.zero();
    vt.zero();
    v.zero();
    
    c = this -> getColumnV(0);
    vt = c + I.getRowV(0)*(c(0)/fabs(c(0)))*(!c);
    v = vt.transpose();
    H = I - (v*vt)*2/(~vt);
    R = H * (* this);
    
    for (int i=1; i<n; i++) {
        R.a[i*m] = (T)0.0;
    }
    Q = H;
    
    //  Other columns
    for (int i=1; i<n-1; i++) {
        c.reSize(n-i);
        vt.reSize(n-i);
        A.reSize(n-i, n-i);
        H.reSize(n-i, n-i);
        v.reSize(n-i, 1);
        
        A = R.getMatrix(i, i);
        c = A.getColumnV(0);
        
        if (c.null()) {
            H.identity();
        } else {
            vt = c + I.getRowV(i, i)*(c(0)/fabs(c(0)))*(!c);
            v = vt.transpose();
            H = I.getMatrix(i, i) - (v*vt)*2/(~vt);
            A = I;
            A.setMatrix(i, i, H);
            H = A;
        }
        
        Q = H*Q;
        R = H*R;
        
        for (int j=0; j<i; j++) {
            R.a[i*m + j] = (T)0.0;
        }
    }
    
    Q = Q.transpose();
    
    for (int i=1; i<n; i++) {
        for (int j=0; j<i; j++) {
            R.a[i*m + j] = (T)0.0;
        }
    }
    
    if ((QorR& Qmatrix) != 0)
        return Q;
    else if ((QorR& Rmatrix) != 0)
        return R;
    else {
        matrixT QR(n,2*n);
        QR.zero();
        QR = Q&&R;
        return QR;
    }
}

MAT_TEMPLATE inline matrixT
matrixT::eigenVectors(int maxIte, unsigned char opt)
{
    if (n != m)
        cout << RM << "eigenVectors()] -> ";
    
    int k;
    T eVa, maxS = 0.0, fact = 1E-07;
    vectorT eVas(n);
    matrixT eVes(n,n), eVe(n,1), Ainv(n,n);
    
    Ainv.zero();
    eVas.zero();
    eVes.zero();
    
    eVas = this -> eigenValues(maxIte, preDef, opt);
    if ((opt& ORGANIZED) != 0)
        eVas = eVas.organize();
    
    if (eVas.duplicate(fact))
        cout << RM << "eigenVectors()] - Hay autovalores duplicados, no se asegura la obtención de todos los autovectores.\n";
    
    for (int i=0; i<n; i++) {
        eVe.ones();
        eVa = (T)eVas(i)+1.0;
        k=1;
        Ainv = (* this) - identity(n,n)*eVas(i)*(1.0+fact);
        Ainv = Ainv^(-1);
        
        while ((eVa < eVas(i)*(1.0-fact) || eVa > eVas(i)*(1.0+fact)) && k < maxIte) {
            eVe = Ainv * eVe;
            maxS = eVe.maxAbs_sig();
            eVa = 1/maxS + eVas(i)*(1.0+fact);
            eVe /= maxS;
            k++;
            
            if ((eVa < eVas(i)*(1.0-fact) || eVa > eVas(i)*(1.0+fact)) && k == maxIte) {
                eVe.ones();
                eVe(i) *= 2.0;
                eVa = (T)eVas(i)+1.0;
                k=1;
                Ainv = (* this) - identity(n,n)*eVas(i)*(1.0+fact);
                Ainv = Ainv^(-1);
                
                while (1/eVa != eVas(i) && k < maxIte) {
                    eVe = Ainv * eVe;
                    maxS = eVe.maxAbs_sig();
                    eVa = 1/maxS + eVas(i)*(1.0+fact);
                    eVe /= maxS;
                    k++;
                }
            }
        }
        eVes.setColumn(i, eVe);
        
        if ((opt& EVe_Ite) != 0)
            cout << RM << "eigenVectors()] - Iteraciones para calcular el autovector " << i+1 << ": " << k << endl;
    }
    return eVes;
}

MAT_TEMPLATE inline matrixT
matrixT::eigenVectors(unsigned char opt)
{
    return (this -> eigenVectors(tolMax, opt));
}

MAT_TEMPLATE inline matrixT
matrixT::eigenVectors()
{
    return (this -> eigenVectors(tolMax, 0));
}

MAT_TEMPLATE inline vectorT
matrixT::eigenVector(T eigenValue, int maxIte, unsigned char opt)
{
    if (n != m)
        cout << RM << "eigenVectors()] -> ";
    
    int k;
    T eVa, maxS = 0.0, fact = 1E-07;
    vectorT eVeV(n);
    matrixT eVe(n,1), Ainv(n,n);
    
    Ainv.zero();
    eVe.ones();
    
    eVa = eigenValue+1.0;
    k=1;
    Ainv = (* this) - identity(n,n)*eigenValue*(1.0+fact);
    Ainv = Ainv^(-1);
    
    while ((eVa < eigenValue*(1.0-fact) || eVa > eigenValue*(1.0+fact)) && k < maxIte) {
        eVe = Ainv * eVe;
        maxS = eVe.maxAbs_sig();
        eVa = 1/maxS + eigenValue*(1.0+fact);
        eVe /= maxS;
        k++;
        
        if ((eVa < eigenValue*(1.0-fact) || eVa > eigenValue*(1.0+fact)) && k == maxIte) {
            eVe.ones();
            eVe(0) *= 2.0;
            eVa = eigenValue+1.0;
            k=1;
            Ainv = (* this) - identity(n,n)*eigenValue*(1.0+fact);
            Ainv = Ainv^(-1);
            
            while (1/eVa != eigenValue && k < maxIte) {
                eVe = Ainv * eVe;
                maxS = eVe.maxAbs_sig();
                eVa = 1/maxS + eigenValue*(1.0+fact);
                eVe /= maxS;
                k++;
            }
        }
    }
    eVeV = eVe.getColumnV(0);
    
    if ((opt& EVe_Ite) != 0)
        cout << RM << "eigenVectors()] - Iteraciones para calcular el autovector: " << k << endl;
    
    return eVeV;
}

MAT_TEMPLATE inline vectorT
matrixT::eigenVector(T eigenValue, unsigned char opt)
{
    return (this -> eigenVector(eigenValue, tolMax, opt));
}

MAT_TEMPLATE inline vectorT
matrixT::eigenVector(T eigenValue)
{
    return (this -> eigenVector(eigenValue, tolMax, 0));
}

MAT_TEMPLATE inline vectorT
matrixT::eigenValues(int maxIte, T factor, unsigned char opt)
{
    if (n != m)
        cout << RM << "eigenValues(maxIte)] -> ";
    
    int k=1;
    T cota = 1.0;
    T err = 0.1;
    matrixT A(n,n), Q(n,2*n), R(n,n);
    vectorT diagOld(n);
    
    A.zero();
    Q.zero();
    R.zero();
    diagOld.zero();
    A = (* this);
    
    while (cota > err && k<maxIte) {
        Q = A.QR(QRmatrix);
        Q||R;
        diagOld = A.getDiagonal();
        A = R*Q;
        cota = (diagOld - A.getDiagonal()).maxAbs();
        err = (A.getDiagonal()).maxAbs()*factor;
        k++;
    }
    
    if ((opt& EVa_Ite) != 0)
        cout << RM << "eigenValues()] - Iteraciones para calcular los autovalores: " << k << endl;
    
    return A.getDiagonal();
}

MAT_TEMPLATE inline vectorT
matrixT::eigenValues(int maxIte, unsigned char opt)
{
    return (this -> eigenValues(maxIte, preDef, opt));
}

MAT_TEMPLATE inline vectorT
matrixT::eigenValues(unsigned char opt)
{
    return (this -> eigenValues(tolMax, preDef, opt));
}

MAT_TEMPLATE inline vectorT
matrixT::eigenValues()
{
    return (this -> eigenValues(tolMax, preDef, 0));
}



//  --- OPERATORS ---
MAT_TEMPLATE inline matrixT&
matrixT::operator=(const T* array)
{
    for (int k=0; k<size; k++) {
        a[k] = array[k];
    }
    return (* this);
}

MAT_TEMPLATE inline matrixT
matrixT::operator+(const matrixT& M)
{
    if (n != M.n || m != M.m)
        cout << RM << "operator+M] - Las matrices no tienen las mismas dimensiones.\n";
    matrixT tmp(n,m);
    tmp.zero();
    for (int k=0; k<size; k++) {
        tmp.a[k] = a[k] + M.a[k];
    }
    return tmp;
}

MAT_TEMPLATE inline matrixT&
matrixT::operator+=(const matrixT& M)
{
    if (n != M.n || m != M.m)
        cout << RM << "operator+=M] - Las matrices no tienen las mismas dimensiones.\n";
    for (int k=0; k<size; k++) {
        a[k] += M.a[k];
    }
    return (* this);
}

MAT_TEMPLATE inline matrixT
matrixT::operator-(const matrixT& M)
{
    if (n != M.n || m != M.m)
        cout << RM << "operator-M] - Las matrices no tienen las mismas dimensiones.\n";
    matrixT tmp(n,m);
    tmp.zero();
    for (int k=0; k<size; k++) {
        tmp.a[k] = a[k] - M.a[k];
    }
    return tmp;
}

MAT_TEMPLATE inline matrixT&
matrixT::operator-=(const matrixT& M)
{
    if (n != M.n || m != M.m)
        cout << RM << "operator-=M] - Las matrices no tienen las mismas dimensiones.\n";
    for (int k=0; k<size; k++) {
        a[k] -= M.a[k];
    }
    return (* this);
}

MAT_TEMPLATE inline matrixT
matrixT::operator*(const T& D)
{
    matrixT tmp(n,m);
    tmp.zero();
    for (int k=0; k<size; k++) {
        tmp.a[k] = a[k] * D;
    }
    return tmp;
}

MAT_TEMPLATE inline matrixT&
matrixT::operator*=(const T& D)
{
    for (int k=0; k<size; k++) {
        a[k] *= D;
    }
    return (* this);
}

MAT_TEMPLATE inline matrixT
matrixT::operator/(const T& D)
{
    matrixT tmp(n,m);
    tmp.zero();
    for (int k=0; k<size; k++) {
        tmp.a[k] = a[k] / D;
    }
    return tmp;
}

MAT_TEMPLATE inline matrixT&
matrixT::operator/=(const T& D)
{
    for (int k=0; k<size; k++) {
        a[k] /= D;
    }
    return (* this);
}

MAT_TEMPLATE inline matrixT
matrixT::operator-()
{
    matrixT tmp(n,m);
    tmp.zero();
    for (int k=0; k<size; k++) {
        tmp.a[k] = -a[k];
    }
    return tmp;
}

MAT_TEMPLATE inline matrixT
matrixT::operator*(const matrixT& M)
{
    matrixT tmp(n,M.m);
    tmp.zero();
    if (m == M.n) {
        for (int i=0; i<n; i++) {
            for (int j=0; j<M.m; j++) {
                for (int k=0; k<m; k++) {
                    tmp.a[i*M.m + j] += a[i*m + k] * M.a[k*M.m + j];
                }
            }
        }
    } else {
        cout << RM << "operator*M] - Las dimensiones de las matrices no permiten su producto.\n";
        tmp.ones();
    }
    return tmp;
}

MAT_TEMPLATE inline matrixT&
matrixT::operator*=(const matrixT& M)
{
    if (n != m || M.n != M.m || n != M.n) {
        cout << RM << "operator*=M] - No se puede multiplicar sobre sí misma porque no es una matriz cuadrada.\n";
        return (* this);
    }
    T prod;
    matrixT tmp(n,M.m);
    tmp.zero();
    for (int i=0; i<n; i++) {
        for (int j=0; j<M.m; j++) {
            prod = (T)0.0;
            for (int k=0; k<M.m; k++) {
                prod += a[i*m + k] * M.a[k*M.m + j];
            }
            tmp.a[i*M.m + j] = prod;
        }
    }
    *this = tmp;
    return (* this);
}

MAT_TEMPLATE inline matrixT
matrixT::operator*(const vectorT& V)
{
    matrixT tmp(n,V.dim());
    tmp.zero();
    if (m != 1) {
        cout << RM << "operator*V] - La matriz debe tener una única columna para multiplicarse por un vector.\n";
        tmp.ones();
    } else {
        for (int i=0; i<n; i++) {
            for (int j=0; j<V.dim(); j++) {
                tmp.a[i*V.dim() + j] = a[i] * V(j);
            }
        }
    }
    return tmp;
}

MAT_TEMPLATE inline matrixT
matrixT::operator^(const int exp)
{
    matrixT tmp(n,m);
    tmp.zero();
    if (n != m) {
        cout << RM << "operator^exp] - No se puede elevar la matriz a " << exp << " porque no es cuadrada.\n";
        tmp.ones();
        return tmp;
    }
    
    if (exp == -1) {
        vectorT I(n), sol(n);
        sol.zero();
        for (int i=0; i<m; i++) {
            I.zero();
            I(i) = (T)1.0;
            sol = this -> solveLU(I);
            for (int k=0; k<m; k++) {
                tmp.a[k*m + i] = sol(k);
            }
            sol.zero();
        }
    }
    
    if (exp > 0) {
        tmp = *this;
        for (int i=2; i<=exp; i++) {
            tmp *= (* this);
        }
    }
    
    return tmp;
}

//  --- MORE OPERATORS ---
MAT_TEMPLATE inline matrixT
operator*(const T& D, const matrixT& M)
{
    matrixT tmp(M.rows(),M.columns());
    tmp.zero();
    for (int k=0; k<M.elements(); k++) {
        tmp(k) = M(k) * D;
    }
    return tmp;
}

MAT_TEMPLATE inline matrixT
operator&&(const matrixT& Mi, const matrixT& Md)
{
    if (Mi.rows() != Md.rows())
        cout << RM << "operator&&(Mi, Md)] - Estas matrices no tienen la misma altura.\n";
    int n,m;
    n = max(Mi.rows(), Md.rows());
    m = Mi.columns() + Md.columns();
    matrixT tmp(n,m);
    tmp.zero();
    
    for (int i=0; i<Mi.rows(); i++) {
        for (int j=0; j<Mi.columns(); j++) {
            tmp(i,j) = Mi(i,j);
        }
    }
    
    for (int i=0; i<Md.rows(); i++) {
        for (int j=Mi.columns(); j<m; j++) {
            tmp(i,j) = Md(i,j-Mi.columns());
        }
    }
    
    return tmp;
}

MAT_TEMPLATE inline void
operator||(matrixT& Mi, matrixT& Md)
{
    if (Mi.rows() != Md.rows())
        cout << RM << "operator||(Mi, Md)] - Las matrices no tienen la misma altura.\n";
    
    Md = Mi.getMatrix(0, (Mi.columns() - Md.columns()), Md.rows(), Md.columns());
    Mi = Mi.getMatrix(0, 0, Mi.rows(), (Mi.columns() - Md.columns()));
}

MAT_TEMPLATE inline void
operator>>(istream& in, matrixT& M)
{
    bool fRow = true;
    unsigned int n=0, m=1, k=0, endR=0;
    string c, cc, num;
    
    while (! in.eof()) {
        getline(in, c);
        cc += c;
        for (int i=0; i<c.size(); i++)
            if (fRow && c[i] == ',')
                m++;
        fRow = false;
        n++;
    }
    
    M.reSize(n, m);
    for (int i=0; i<n-1; i++) {
        endR = cc.find("\r", endR+1);
        for (int j=0; j<m; j++) {
            num = "";
            while (cc[k] != ',' && k < endR) {
                num += cc[k];
                k++;
            }
            k++;
            M(i,j) = strtod(num.c_str(), NULL);
        }
        k = endR + 1;
    }
    
    for (int j=0; j<m; j++) {
        num = "";
        while (cc[k] != ',' && (cc[k] != '\r' && cc[k] != '\n')) {
            num += cc[k];
            k++;
            if (k == cc.size())
                break;
        }
        k++;
        M(n-1,j) = strtod(num.c_str(), NULL);
    }
}

MAT_TEMPLATE inline ostream&
operator<<(ostream& out, const matrixT& M)
{
//    if (out == cout) {
//        out.width();
//        out << fixed;
//        out.fill(' ');
//        out.precision(6);
//        if (M.rows() == 1 && M.columns() == 1) {
//            out << M(0) << endl;
//        } else if (M.rows() == 1) {
//            out << "[" << M(0);
//            for (int j=1; j<M.columns(); j++) {
//                out << "\t" << M(j); 
//            }
//            out << "]\n";
//        } else {
//            out << "⎡" << M(0);
//            for (int j=1; j<M.columns(); j++) {
//                out << right << "\t" << M(j);
//            }
//            out.width();
//            out << "\t⎤\n";
//            for (int i=1; i<M.rows()-1; i++) {
//                out << "⎢" << M(i*M.columns());
//                for (int j=1; j<M.columns(); j++) {
//                    out << "\t" << M(i*M.columns() + j);
//                }
//                out << "\t⎥\n";
//            }
//            out << "⎣" << M((M.rows()-1)*M.columns());
//            for (int j=1; j<M.columns(); j++) {
//                out << "\t" << M((M.rows()-1)*M.columns() + j);
//            }
//            out << "\t⎦\n";
//        }
//        out.precision();
//    } else {
//        for (int i=0; i<M.rows(); i++) {
//            for (int j=0; j<M.columns()-1; j++) {
//                out << M(i*M.columns() + j) << ",";
//            }
//            out << M(i*M.columns() + (M.columns()-1));
//            if (i != M.rows()-1)
//                out << endl; 
//        }
//    }
    return out;
}



//  --- OTHER FUNCTIONS ---
MAT_TEMPLATE inline matrixT
zero(int _n, int _m)
{
    matrixT tmp(_n,_m);
    tmp.zero();
    return tmp;
}

MAT_TEMPLATE inline matrixT
ones(int _n, int _m)
{
    matrixT tmp(_n,_m);
    tmp.ones();
    return tmp;
}

MAT_TEMPLATE inline matrixT
identity(int _n, int _m)
{
    if (_n != _m)
        cout << RM << "identity(n, m)] - No es una matriz cuadrada.\n";
    matrixT I(_n,_m);
    I.identity();
    return I;
}

MAT_TEMPLATE inline matrixT
setdiff(matrixT& A, const matrixT& B, const int& reps)
{
    matrixT tmp(A.rows(), A.columns());
    int dim = 0, coincidence = 0;
    
    for (int i=0; i<A.rows(); i++) {
        
        bool equal = false;
        for (int j=0; j<B.rows(); j++) {
            
            int count=0;
            for (int k=0; k<A.columns(); k++) {
                if (A(i,k) == B(j,k))
                    count++;
            }
            
            if (count == A.columns()) {
                equal = true;
                coincidence++;
                break;
            }
        }
        
        if (equal != true) {
            tmp.setRowV(dim,A.getRowV(i));
            dim++;
        }
        
        if (coincidence == reps) {
            tmp.setMatrix(dim,0,A.getMatrix(i+1,0));
            dim = A.rows() - reps;
            break;
        }
    }
    
    if (dim == 0) {
        matrixT fail(1,1);
        fail(0,0) = 0.0;
        return fail;
    } else {
        tmp.reSize(dim, A.columns());
        return tmp;
    }
}

MAT_TEMPLATE inline matrixT
setdiff(matrixT& A, const matrixT& B)
{
    return setdiff(A, B, A.rows());
}
