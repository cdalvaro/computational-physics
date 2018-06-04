//
//  matrix.cpp
//  Class for vectorial operations
//  Vibration Normal Modes 2D
//
//  Created by Carlos David on 5/8/12.
//  Copyright (c) 2012 NEPTUNO. All rights reserved.
//

#include "matrix.hpp"

#include <typeinfo>


const int tolMax = 1000;
const double preDef = 1E-04;
const std::string RM = "\n[Matrix::";

using namespace cda::math::containers;

            //  --- MATRIX CLASS ---

//  --- DEFINITION ---

template<typename T>
Matrix<T>::Matrix(T* _a)
{
    a = _a;
}

template<typename T>
Matrix<T>::Matrix(const Matrix<T>& M)
{
    n = M.n;
    m = M.m;
    size = M.size;
    a = new T[size];
    std::copy(M.a, M.a + M.size, a);
}

template<typename T>
template <class T2>
Matrix<T>::Matrix(const Matrix<T2>& M)
{
    n = M.rows();
    m = M.columns();
    size = M.elements();
    a = new T[size];
    
    for (int k=0; k<size; k++) {
        a[k] = (T)M(k);
    }
}

template<typename T>
Matrix<T>::~Matrix()
{
    delete a;
    a = nullptr;
}

template<typename T> inline void
Matrix<T>::reSize(int _n, int _m)
{
    Matrix<T> tmp(n,m);
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

template<typename T> inline
Matrix<T> &Matrix<T>::operator=(const Matrix<T>& M)
{
    if (&M != this) {
        delete a;
        n = M.n;
        m = M.m;
        size = M.size;
        a = new T[size];
        std::copy(M.a, M.a + M.size, a);
    }
    return (* this);
}



//  --- GETS AND SETS ---
//  Gets
template<typename T> inline Matrix<T>
Matrix<T>::getRow(int _i)
{
    Matrix<T> tmp(1,m);
    tmp.zero();
    if (_i >= n) {
        std::cout << RM << "getRow(i)] - Esta fila no existe, se ha devuelvo una fila de unos.\n";
        tmp.ones();
    } else {
        for (int j=0; j<m; j++) {
            tmp.a[j] = a[_i*m + j];
        }
    }
    return tmp;
}

template<typename T> inline Matrix<T>
Matrix<T>::getColumn(int _j)
{
    Matrix<T> tmp(n,1);
    tmp.zero();
    if (_j >= m) {
        std::cout << RM << "getColumn(j)] - Esta columna no existe, se ha devuelvo una columna de unos.\n";
        tmp.ones();
    } else {
        for (int i=0; i<n; i++) {
            tmp.a[i] = a[i*m + _j];
        }
    }
    return tmp;
}

template<typename T> inline Matrix<T>
Matrix<T>::getMatrix(int _i, int _j, int height, int columns)
{
    Matrix<T> tmp(height,columns);
    tmp.zero();
    if (n < _i+height || m < _j+columns) {
        std::cout << RM << "getMatrix(i, j, height, columns)] - Estás intentado acceder a una parte de la matriz que no existe, se ha devuelto una matriz de unos\n";
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

template<typename T> inline Matrix<T>
Matrix<T>::getMatrix(int _i, int _j)
{
    return (this -> getMatrix(_i, _j, n-_i, m-_j));
}

template<typename T> inline Vector<T>
Matrix<T>::getDiagonal()
{
    Vector<T> tmp(n);
    if (n != m)
        std::cout << RM << "getDiagonal()] - Esta no es una matriz cuadrada.\n";
    for (int i=0; i<n; i++) {
        tmp[i] = a[i*m + i];
    }
    return tmp;
}

template<typename T> inline Vector<T>
Matrix<T>::getRowV(int _i, int _j)
{
    Vector<T> tmp(m-_j);
    if (_i >= n) {
        std::cout << RM << "getRowV(i, j)] - Esta fila no existe, se ha devuelvo una fila de unos.\n";
        tmp.Ones();
    } else {
        for (int j=0; j<m-_j; j++) {
            tmp[j] = a[_i*m + (j+_j)];
        }
    }
    return tmp;
}

template<typename T> inline Vector<T>
Matrix<T>::getRowV(int _i)
{
    return (this -> getRowV(_i, 0));
}

template<typename T> inline Vector<T>
Matrix<T>::getColumnV(int _i, int _j)
{
    Vector<T> tmp(n-_i);
    if (_j >= m) {
        std::cout << RM << "getColumnV(i, j)] - Esta columna no existe, se ha devuelvo una columna de unos.\n";
        tmp.Ones();
    } else {
        for (int i=0; i<n-_i; i++) {
            tmp[i] = a[(i+_i)*m + _j];
        }
    }
    return tmp;
}

template<typename T> inline Vector<T>
Matrix<T>::getColumnV(int _j)
{
    return (this -> getColumnV(0, _j));
}


//  Sets
template<typename T> inline void
Matrix<T>::setColumn(int _i, int _j, const Matrix<T>& Mc)
{
    if (n < _i+Mc.n || Mc.m > 1)
        std::cout << RM << "setColumn(i, j, Mc)] - La columna que se está intentando sustituir es demasiado larga o no se está introduciendo una matriz columna.\n";
    for (int i=_i; i<_i+Mc.n; i++) {
        a[i*m + _j] = Mc.a[(i-_i)*Mc.m];
    }
}

template<typename T> inline void
Matrix<T>::setColumn(int _j, const Matrix<T>& Mc)
{
    this -> setColumn(0, _j, Mc);
}

template<typename T> inline void
Matrix<T>::setColumnV(int _i, int _j, const Vector<T>& V)
{
    if (n < _i+V.Size())
        std::cout << RM << "setColumnV(i, j, V)] - La columna que se está intentando sustituir es demasiado larga.\n";
    for (int i=_i; i<_i+V.Size(); i++) {
        a[i*m + _j] = V[i-_i];
    }
}

template<typename T> inline void
Matrix<T>::setColumnV(int _j, const Vector<T>& V)
{
    this -> setColumnV(0, _j, V);
}

template<typename T> inline void
Matrix<T>::setRow(int _i, int _j, const Matrix<T>& Mr)
{
    if (m < _j+Mr.m || Mr.n > 1) {
        std::cout << RM << "setRow(i, j, Mr)] - La fila que se está intentando sustituir es demasiado larga o no se está introduciendo una matriz fila.\n";
    }
    for (int j=_j; j<_j+Mr.m; j++) {
        a[_i*m + j] = Mr.a[j-_j];
    }
}

template<typename T> inline void
Matrix<T>::setRow(int _i, const Matrix<T>& Mr)
{
    this -> setRow(_i, 0, Mr);
}

template<typename T> inline void
Matrix<T>::setRowV(int _i, int _j, const Vector<T>& V)
{
    if (m < _j+V.Size())
        std::cout << RM << "setRowV(i, j, V)] - La fila que se está intentando sustituir es demasiado larga.\n";
    for (int j=_j; j<_j+V.Size(); j++) {
        a[_i*m + j] = V[j-_j];
    }
}

template<typename T> inline void
Matrix<T>::setRowV(int _i, const Vector<T>& V)
{
    this -> setRowV(_i, 0, V);
}

template<typename T> inline void
Matrix<T>::setMatrix(int _i, int _j, const Matrix<T>& M)
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
template<typename T> inline int
Matrix<T>::dims() const
{
    int dims[2];
    dims[0] = n;
    dims[1] = m;
    return (* dims);
}

template<typename T> inline int
Matrix<T>::rows() const
{
    return n;
}

template<typename T> inline int
Matrix<T>::columns() const
{
    return m;
}

template<typename T> inline int
Matrix<T>::elements() const
{
    return size;
}


template<typename T> inline T
Matrix<T>::sumRow(int _i) const
{
    if (n >= _i)
        std::cout << RM << "sumRow(i)] - Esta fila no existe.\n";
    T sum = (T)0.0;
    for (int j=0; j<m; j++) {
        sum += a[_i*m + j];
    }
    return sum;
}

template<typename T> inline T
Matrix<T>::sumColumn(int _j) const
{
    if (m >= _j)
        std::cout << RM << "sumColumn(j)] - Esta columna no existe.\n";
    T sum = (T)0.0;
    for (int i=0; i<n; i++) {
        sum += a[i*m + _j];
    }
    return sum;
}

template<typename T> inline T
Matrix<T>::max() const
{
    T max = a[0];
    for (int k=1; k<size; k++) {
        if (max < a[k])
            max = a[k];
    }
    return max;
}

template<typename T> inline T
Matrix<T>::maxAbs() const
{
    T max = fabs(a[0]);
    for (int k=1; k<size; k++) {
        if (fabs(a[k]) > max)
            max = fabs(a[k]);
    }
    return max;
}

template<typename T> inline T
Matrix<T>::maxAbs_sig() const
{
    T max = a[0];
    for (int k=1; k<size; k++) {
        if (fabs(a[k]) > fabs(max))
            max = a[k];
    }
    return max;
}

template<typename T> inline T
Matrix<T>::min() const
{
    T min = a[0];
    for (int k=1; k<size; k++) {
        if (min > a[k])
            min = a[k];
    }
    return min;
}

template<typename T> inline T
Matrix<T>::det()
{
    T det = 1.0;
    if (n != m) {
        std::cout << RM << "det()] - La matriz no es cuadrada, así que no se puede calcular su determinante.\n";
        return det;
    }
    Matrix<T> tmp(n,m);
    tmp.zero();
    tmp = this -> U();
    
    for (int i=0; i<n; i++) {
        det *= tmp(i,i);
    }
    
    return det;
}

template<typename T> inline T&
Matrix<T>::operator()(int _i, int _j) const
{
    if (n <= _i || m <= _j)
        std::cout << RM << "operator(i,j)] - Este elemento de matriz no existe.\n";
    return a[_i*m + _j];
}

template<typename T> inline T&
Matrix<T>::operator()(int _k) const
{
    if (size <= _k)
        std::cout << RM << "operator(k)] - Este elemento de matriz no existe.\n";
    return a[_k];
}


//  Returns a matrix
template<typename T> inline Matrix<T>
Matrix<T>::sumRows()
{
    T sum;
    Matrix<T> tmp(n,1);
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

template<typename T> inline Matrix<T>
Matrix<T>::sumColumns()
{
    T sum;
    Matrix<T> tmp(1,m);
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

template<typename T> inline Matrix<T>
Matrix<T>::transpose()
{
    Matrix<T> tmp(m,n);
    tmp.zero();
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            tmp.a[i*n + j] = a[j*m + i];
        }
    }
    return tmp;
}

template<typename T> inline Matrix<T>
Matrix<T>::zero(int _n, int _m)
{
    Matrix<T> tmp(_n,_m);
    tmp.zero();
    return tmp;
}

template<typename T> inline Matrix<T>
Matrix<T>::ones(int _n, int _m)
{
    Matrix<T> tmp(_n,_m);
    tmp.ones();
    return tmp;
}

template<typename T> inline Matrix<T>
Matrix<T>::identity(int _n, int _m)
{
    if (_n != _m)
        std::cout << RM << "identity(n, m)] - No es una matriz cuadrada.\n";
    Matrix<T> I(_n,_m);
    I.identity();
    return I;
}


template<typename T> inline Vector<T>
Matrix<T>::sumRowsV()
{
    T sum;
    Vector<T> tmp(n);
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

template<typename T> inline Vector<T>
Matrix<T>::sumColumnsV()
{
    T sum;
    Vector<T> tmp(m);
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
template<typename T> inline bool
Matrix<T>::null()
{
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            if (a[i*m + j] != (T)0.0)
                return false;
        }
    }
    return true;
}

template<typename T> inline bool
Matrix<T>::duplicate(T range)
{
    for (int k=0; k<size; k++) {
        for (int p=0; p<size; p++) {
            if (k != p && (a[p] >= a[k]*(1.0-range) && a[p] <= a[k]*(1.0+range)))
                return true;
        }
    }
    return false;
}

template<typename T> inline bool
Matrix<T>::duplicate()
{
    return (this -> duplicate(0.0));
}


//  Void functions
template<typename T> inline void
Matrix<T>::zero()
{
    for (int k=0; k<size; k++) {
        if (typeid(T) == typeid(bool))
            a[k] = false;
        else
            a[k] = (T)0.0;
    }
}

template<typename T> inline void
Matrix<T>::ones()
{
    for (int k=0; k<size; k++) {
        if (typeid(T) == typeid(bool))
            a[k] = true;
        else
            a[k] = (T)1.0;
    }
}

template<typename T> inline void
Matrix<T>::identity()
{
    if (n != m)
        std::cout << RM << "identity()] - No es una matriz cuadrada.\n";
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

template<typename T> inline void
Matrix<T>::write(const std::string& path, const std::string& filename)
{
    const std::string Npath = path + filename;
    std::ofstream out(Npath.data());
    out << (* this);
}

template<typename T> inline void
Matrix<T>::write(const std::string& filename)
{
    const std::string home = getenv("HOME");
    const std::string path = home + "/Desktop/";
    write(path, filename);
}



//  --- MÉTODO LU ---
template<typename T> inline Matrix<T>
Matrix<T>::LU()
{
    if (n != m) {
        std::cout << RM << "LU()] - No se puede calcular la diagonal inferior porque no es una matriz cuadrada.\n";
        Matrix<T> tmp(n,m);
        tmp.ones();
        return tmp;
    }
    
    Matrix<T> LU(n,m);
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

template<typename T> inline Matrix<T>
Matrix<T>::U()
{
    if (n != m)
        std::cout << RM << "U()] -> ";
    Matrix<T> LU(n,m), U(n,m);
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

template<typename T> inline Matrix<T>
Matrix<T>::L()
{
    if (n != m)
        std::cout << RM << "L()] -> ";
    Matrix<T> LU(n,m), L(n,m);
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

template<typename T> inline Vector<T>
Matrix<T>::solveLU(const Vector<T>& B)
{
    Vector<T> tmp(n, 0);
    if (n != m || n != B.Size()) {
        std::cout << RM << "solveLU(B)] - El sistema no es compatible determinado, la matriz no es cuadrada o el vector de términos independientes no coincide con el número de incógnitas.\n";
        tmp.Ones();
        return tmp;
    }
    Vector<T> x(n, 0);
    Matrix<T> LU(n,m);
    LU.zero();
    LU = this -> LU();
    
    T sum = (T)0.0;
    for (int i=0; i<n; i++) {
        sum = (T)0.0;
        for (int j=0; j<i; j++) {
            sum += LU.a[i*m + j] * tmp[j];
        }
        tmp[i] = B[i] - sum;
    }
    
    for (int i=n-1; i>=0; i--) {
        sum = (T)0.0;
        for (int j=m-1; j>=0; j--) {
            sum += LU.a[i*m + j] * x[j];
        }
        x[i] = (tmp[i]-sum) / LU.a[i*m + i];
    }
    return x;
}

template<typename T> inline Vector<T>
Matrix<T>::solveLU3d(const Vector<T>& B)
{
    Vector<T> tmp(n, 0);
    if (n != B.Size()) {
        std::cout << RM << "solveLU3d(B)] - El sistema no es compatible determinado, la matriz no es cuadrada o el vector de términos independientes no coincide con el número de incógnitas.\n";
        tmp.Ones();
        return tmp;
    }
    Vector<T> al(n, 0), be(n, 0), x(n, 0);
    
    be[0] = a[1];
    for (int i=1; i<n; i++) {
        al[i] = a[i*m] / be[i-1];
        be[i] = a[i*m + 1] - al[i]*a[(i-1)*m + 2];
    }
    tmp[0] = B[0];
    for (int i=1; i<n; i++) {
        tmp[i] = B[i] - al[i]*tmp[i-1];
    }
    x[n-1] = tmp[n-1]/be[n-1];
    for (int i=(n-2); i>=0; i--) {
        x[i] = (tmp[i] - a[i*m + 2]*x[i+1]) / be[i];
    }
    
    return x;
}

template<typename T> inline Vector<T>
Matrix<T>::solveGS3d(const Vector<T>& B, T err)
{
    Vector<T> x(n, 0);
    if (n != B.Size()) {
        std::cout << RM << "solveGS3d(B, err)] - El sistema no es compatible determinado, la matriz no es cuadrada o el vector de términos independientes no coincide con el número de incógnitas.\n";
        x.Ones();
        return x;
    }
    
    T error = 0.0, xAnt;
    do {
        error = 0.0;
        
        xAnt = x[0];
        x[0] = (B[0] - a[0]*x[1])/a[1];
        error += std::pow((x[0]-xAnt), 2);
        
        for (int i=1; i<(n-1); i++) {
            xAnt = x[i];
            x[i] = (B[i] - a[i*n]*x[i-1] - a[i*n + 2]*x[i+1])/a[i*n + 1];
            error += std::pow((x[i]-xAnt), 2);
        }
        
        xAnt = x[n-1];
        x[n-1] = (B[n-1] - a[(n-1)*n]*x[n-2])/a[(n-1)*n + 1];
        error += std::pow((x[n-1]-xAnt), 2);
        
        error = std::sqrt(error);
        
    } while (error > err);
    return x;
}


//  --- MÉTODO QR ---
template<typename T> inline Matrix<T>
Matrix<T>::QR(unsigned char QorR)
{
    if (n != m) {
        std::cout << RM << "QR(QorR) - La matriz no es cuadrada.";
        return ones(n, n);
    }
    
    Matrix<T> Q(n,n), R(n,n), I(n,n), H(n,n), A(n,n);
    Q.zero();
    R.zero();
    H.zero();
    A.zero();
    I.identity();
    
    //  First column
    Vector<T> c(n, 0), vt(n, 0);
    Matrix<T> v(n,1);
    v.zero();
    
    c = this->getColumnV(0);
    vt = c + I.getRowV(0)*(c[0]/std::abs(c[0]))*(c.Norm());
    // TODO: Implement transpose method for Vector class
    // v = vt.transpose();
    for (auto it = vt.Begin(); it != vt.End(); ++it) {
        v(std::distance(vt.Begin(), it), 0) = *it;
    }
    
    H = I - (v*vt)*2/(vt.SquaredNorm());
    R = H * (* this);
    
    for (int i=1; i<n; i++) {
        R.a[i*m] = (T)0.0;
    }
    Q = H;
    
    //  Other columns
    for (int i=1; i<n-1; i++) {
        c.Resize(n-i);
        vt.Resize(n-i);
        A.reSize(n-i, n-i);
        H.reSize(n-i, n-i);
        v.reSize(n-i, 1);
        
        A = R.getMatrix(i, i);
        c = A.getColumnV(0);
        
        if (c.IsNull()) {
            H.identity();
        } else {
            vt = c + I.getRowV(i, i)*(c[0]/std::abs(c[0]))*(c.Norm());
            // TODO: Implement transpose method for Vector class
            // v = vt.transpose();
            for (auto it = vt.Begin(); it != vt.End(); ++it) {
                v(std::distance(vt.Begin(), it), 0) = *it;
            }
            
            H = I.getMatrix(i, i) - (v*vt)*2/(vt.SquaredNorm());
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
        Matrix<T> QR(n,2*n);
        QR.zero();
        QR = Q&&R;
        return QR;
    }
}

template<typename T> inline Matrix<T>
Matrix<T>::eigenVectors(int maxIte, unsigned char opt)
{
    if (n != m)
        std::cout << RM << "eigenVectors()] -> ";
    
    int k;
    T eVa, maxS = 0.0, fact = 1E-07;
    Vector<T> eVas(n);
    Matrix<T> eVes(n,n), eVe(n,1), Ainv(n,n);
    
    Ainv.zero();
    eVas.zero();
    eVes.zero();
    
    eVas = this -> eigenValues(maxIte, preDef, opt);
    if ((opt& ORGANIZED) != 0)
        eVas = eVas.organize();
    
    if (eVas.duplicate(fact))
        std::cout << RM << "eigenVectors()] - Hay autovalores duplicados, no se asegura la obtención de todos los autovectores.\n";
    
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
            std::cout << RM << "eigenVectors()] - Iteraciones para calcular el autovector " << i+1 << ": " << k << std::endl;
    }
    return eVes;
}

template<typename T> inline Matrix<T>
Matrix<T>::eigenVectors(unsigned char opt)
{
    return (this -> eigenVectors(tolMax, opt));
}

template<typename T> inline Matrix<T>
Matrix<T>::eigenVectors()
{
    return (this -> eigenVectors(tolMax, 0));
}

template<typename T> inline Vector<T>
Matrix<T>::eigenVector(T eigenValue, int maxIte, unsigned char opt)
{
    if (n != m)
        std::cout << RM << "eigenVectors()] -> ";
    
    int k;
    T eVa, maxS = 0.0, fact = 1E-07;
    Vector<T> eVeV(n);
    Matrix<T> eVe(n,1), Ainv(n,n);
    
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
        std::cout << RM << "eigenVectors()] - Iteraciones para calcular el autovector: " << k << std::endl;
    
    return eVeV;
}

template<typename T> inline Vector<T>
Matrix<T>::eigenVector(T eigenValue, unsigned char opt)
{
    return (this -> eigenVector(eigenValue, tolMax, opt));
}

template<typename T> inline Vector<T>
Matrix<T>::eigenVector(T eigenValue)
{
    return (this -> eigenVector(eigenValue, tolMax, 0));
}

template<typename T> inline Vector<T>
Matrix<T>::eigenValues(int maxIte, T factor, unsigned char opt)
{
    if (n != m)
        std::cout << RM << "eigenValues(maxIte)] -> ";
    
    int k=1;
    T cota = 1.0;
    T err = 0.1;
    Matrix<T> A(n,n), Q(n,2*n), R(n,n);
    Vector<T> diagOld(n, 0);
    
    A.zero();
    Q.zero();
    R.zero();
    A = (* this);
    
    while (cota > err && k<maxIte) {
        Q = A.QR(QRmatrix);
        Q||R;
        diagOld = A.getDiagonal();
        A = R*Q;
        cota = (diagOld - A.getDiagonal()).AbsoluteMaximumElement();
        err = (A.getDiagonal()).AbsoluteMaximumElement()*factor;
        k++;
    }
    
    if ((opt& EVa_Ite) != 0)
        std::cout << RM << "eigenValues()] - Iteraciones para calcular los autovalores: " << k << std::endl;
    
    return A.getDiagonal();
}

template<typename T> inline Vector<T>
Matrix<T>::eigenValues(int maxIte, unsigned char opt)
{
    return (this -> eigenValues(maxIte, preDef, opt));
}

template<typename T> inline Vector<T>
Matrix<T>::eigenValues(unsigned char opt)
{
    return (this -> eigenValues(tolMax, preDef, opt));
}

template<typename T> inline Vector<T>
Matrix<T>::eigenValues()
{
    return (this -> eigenValues(tolMax, preDef, 0));
}



//  --- OPERATORS ---
template<typename T> inline Matrix<T>&
Matrix<T>::operator=(const T* array)
{
    for (int k=0; k<size; k++) {
        a[k] = array[k];
    }
    return (* this);
}

template<typename T> inline Matrix<T>
Matrix<T>::operator+(const Matrix<T>& M)
{
    if (n != M.n || m != M.m)
        std::cout << RM << "operator+M] - Las matrices no tienen las mismas dimensiones.\n";
    Matrix<T> tmp(n,m);
    tmp.zero();
    for (int k=0; k<size; k++) {
        tmp.a[k] = a[k] + M.a[k];
    }
    return tmp;
}

template<typename T> inline Matrix<T>&
Matrix<T>::operator+=(const Matrix<T>& M)
{
    if (n != M.n || m != M.m)
        std::cout << RM << "operator+=M] - Las matrices no tienen las mismas dimensiones.\n";
    for (int k=0; k<size; k++) {
        a[k] += M.a[k];
    }
    return (* this);
}

template<typename T> inline Matrix<T>
Matrix<T>::operator-(const Matrix<T>& M)
{
    if (n != M.n || m != M.m)
        std::cout << RM << "operator-M] - Las matrices no tienen las mismas dimensiones.\n";
    Matrix<T> tmp(n,m);
    tmp.zero();
    for (int k=0; k<size; k++) {
        tmp.a[k] = a[k] - M.a[k];
    }
    return tmp;
}

template<typename T> inline Matrix<T>&
Matrix<T>::operator-=(const Matrix<T>& M)
{
    if (n != M.n || m != M.m)
        std::cout << RM << "operator-=M] - Las matrices no tienen las mismas dimensiones.\n";
    for (int k=0; k<size; k++) {
        a[k] -= M.a[k];
    }
    return (* this);
}

template<typename T> inline Matrix<T>
Matrix<T>::operator*(const T& D)
{
    Matrix<T> tmp(n,m);
    tmp.zero();
    for (int k=0; k<size; k++) {
        tmp.a[k] = a[k] * D;
    }
    return tmp;
}

template<typename T> inline Matrix<T>&
Matrix<T>::operator*=(const T& D)
{
    for (int k=0; k<size; k++) {
        a[k] *= D;
    }
    return (* this);
}

template<typename T> inline Matrix<T>
Matrix<T>::operator/(const T& D)
{
    Matrix<T> tmp(n,m);
    tmp.zero();
    for (int k=0; k<size; k++) {
        tmp.a[k] = a[k] / D;
    }
    return tmp;
}

template<typename T> inline Matrix<T>&
Matrix<T>::operator/=(const T& D)
{
    for (int k=0; k<size; k++) {
        a[k] /= D;
    }
    return (* this);
}

template<typename T> inline Matrix<T>
Matrix<T>::operator-()
{
    Matrix<T> tmp(n,m);
    tmp.zero();
    for (int k=0; k<size; k++) {
        tmp.a[k] = -a[k];
    }
    return tmp;
}

template<typename T> inline Matrix<T>
Matrix<T>::operator*(const Matrix<T>& M)
{
    Matrix<T> tmp(n,M.m);
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
        std::cout << RM << "operator*M] - Las dimensiones de las matrices no permiten su producto.\n";
        tmp.ones();
    }
    return tmp;
}

template<typename T> inline Matrix<T>&
Matrix<T>::operator*=(const Matrix<T>& M)
{
    if (n != m || M.n != M.m || n != M.n) {
        std::cout << RM << "operator*=M] - No se puede multiplicar sobre sí misma porque no es una matriz cuadrada.\n";
        return (* this);
    }
    T prod;
    Matrix<T> tmp(n,M.m);
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

template<typename T> inline Matrix<T>
Matrix<T>::operator*(const Vector<T>& V)
{
    Matrix<T> tmp(n,V.Size());
    tmp.zero();
    if (m != 1) {
        std::cout << RM << "operator*V] - La matriz debe tener una única columna para multiplicarse por un vector.\n";
        tmp.ones();
    } else {
        for (int i=0; i<n; i++) {
            for (int j=0; j<V.Size(); j++) {
                tmp.a[i*V.Size() + j] = a[i] * V[j];
            }
        }
    }
    return tmp;
}

template<typename T> inline Matrix<T>
Matrix<T>::operator^(const int exp)
{
    Matrix<T> tmp(n,m);
    tmp.zero();
    if (n != m) {
        std::cout << RM << "operator^exp] - No se puede elevar la matriz a " << exp << " porque no es cuadrada.\n";
        tmp.ones();
        return tmp;
    }
    
    if (exp == -1) {
        Vector<T> I(n);
        for (int i=0; i<m; i++) {
            I.Zero();
            I[i] = (T)1.0;
            auto sol(this->solveLU(I));
            for (int k=0; k<m; k++) {
                tmp.a[k*m + i] = sol[k];
            }
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
template<typename T> inline Matrix<T>
operator*(const T& D, const Matrix<T>& M)
{
    Matrix<T> tmp(M.rows(),M.columns());
    tmp.zero();
    for (int k=0; k<M.elements(); k++) {
        tmp(k) = M(k) * D;
    }
    return tmp;
}

template<typename T> inline Matrix<T>
operator&&(const Matrix<T>& Mi, const Matrix<T>& Md)
{
    if (Mi.rows() != Md.rows())
        std::cout << RM << "operator&&(Mi, Md)] - Estas matrices no tienen la misma altura.\n";
    int n,m;
    n = std::max(Mi.rows(), Md.rows());
    m = Mi.columns() + Md.columns();
    Matrix<T> tmp(n,m);
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

template<typename T> inline void
operator||(Matrix<T>& Mi, Matrix<T>& Md)
{
    if (Mi.rows() != Md.rows())
        std::cout << RM << "operator||(Mi, Md)] - Las matrices no tienen la misma altura.\n";
    
    Md = Mi.getMatrix(0, (Mi.columns() - Md.columns()), Md.rows(), Md.columns());
    Mi = Mi.getMatrix(0, 0, Mi.rows(), (Mi.columns() - Md.columns()));
}

template<typename T> inline void
operator>>(std::istream& in, Matrix<T>& M)
{
    bool fRow = true;
    unsigned int n=0, m=1, k=0, endR=0;
    std::string c, cc, num;
    
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

template<typename T> inline std::ostream&
operator<<(std::ostream& out, const Matrix<T>& M)
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
template<typename T>
Matrix<T> cda::math::containers::zero(int _n, int _m) {
    Matrix<T> tmp(_n,_m);
    tmp.zero();
    return tmp;
}

template<typename T> inline Matrix<T>
cda::math::containers::ones(int _n, int _m)
{
    Matrix<T> tmp(_n,_m);
    tmp.ones();
    return tmp;
}

template<typename T> inline Matrix<T>
cda::math::containers::identity(int _n, int _m)
{
    if (_n != _m)
        std::cout << RM << "identity(n, m)] - No es una matriz cuadrada.\n";
    Matrix<T> I(_n,_m);
    I.identity();
    return I;
}

template<typename T> inline Matrix<T>
cda::math::containers::setdiff(Matrix<T>& A, const Matrix<T>& B, const int& reps)
{
    Matrix<T> tmp(A.rows(), A.columns());
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
        Matrix<T> fail(1,1);
        fail(0,0) = 0.0;
        return fail;
    } else {
        tmp.reSize(dim, A.columns());
        return tmp;
    }
}

template<typename T> inline Matrix<T>
cda::math::containers::setdiff(Matrix<T>& A, const Matrix<T>& B)
{
    return setdiff(A, B, A.rows());
}
