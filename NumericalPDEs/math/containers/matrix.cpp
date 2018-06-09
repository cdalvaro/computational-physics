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
    tmp.Zero();
    tmp = this -> U();
    
    for (int i=0; i<n; i++) {
        det *= tmp(i,i);
    }
    
    return det;
}

//  Returns a matrix
template<typename T> inline Matrix<T>
Matrix<T>::sumRows()
{
    T sum;
    Matrix<T> tmp(n,1);
    tmp.Zero();
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
    tmp.Zero();
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
    tmp.Zero();
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            tmp.a[i*n + j] = a[j*m + i];
        }
    }
    return tmp;
}



template<typename T> inline Vector<T>
Matrix<T>::sumRowsV()
{
    T sum;
    Vector<T> tmp(n);
    tmp.Zero();
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
    tmp.Zero();
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
        tmp.Ones();
        return tmp;
    }
    
    Matrix<T> LU(n,m);
    LU.Zero();
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
    LU.Zero();
    U.Zero();
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
    LU.Zero();
    L.Zero();
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
    LU.Zero();
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
        return Ones(n, n);
    }
    
    Matrix<T> Q(n,n), R(n,n), I(n,n), H(n,n), A(n,n);
    Q.Zero();
    R.Zero();
    H.Zero();
    A.Zero();
    I.Identity();
    
    //  First column
    auto c = this->GetColumnAsVector(0);
    auto vt = c + I.GetRowAsVector(0)*(c[0]/std::abs(c[0]))*(c.Norm());
    
    // TODO: Implement transpose method for Vector class
    // v = vt.transpose();
    Matrix<T> v(n, 1);
    for (auto it = vt.Begin(); it != vt.End(); ++it) {
        v[std::distance(vt.Begin(), it)][0] = *it;
    }
    
    H = I - (v*vt) * 2.0/(vt.SquaredNorm());
    R = H * (* this);
    
    for (int i=1; i<n; i++) {
        R.a[i*m] = (T)0.0;
    }
    Q = H;
    
    //  Other columns
    for (int i=1; i<n-1; i++) {
        c.Resize(n-i);
        vt.Resize(n-i);
        A.Resize(n-i, n-i);
        H.Resize(n-i, n-i);
        v.Resize(n-i, 1);
        
        A = R.GetMatrix(i, i);
        c = A.GetColumnAsVector(0);
        
        if (c.IsNull()) {
            H.Identity();
        } else {
            vt = c + I.GetRowAsVector(i, i)*(c[0]/std::abs(c[0]))*(c.Norm());
            // TODO: Implement transpose method for Vector class
            // v = vt.transpose();
            for (auto it = vt.Begin(); it != vt.End(); ++it) {
                v[std::distance(vt.Begin(), it)][0] = *it;
            }
            
            H = I.GetMatrix(i, i) - 2.0*(v*vt)/(vt.SquaredNorm());
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
        QR.Zero();
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
    
    Ainv.Zero();
    eVas.Zero();
    eVes.Zero();
    
    eVas = this -> eigenValues(maxIte, preDef, opt);
    if ((opt& ORGANIZED) != 0)
        eVas = eVas.organize();
    
    if (eVas.duplicate(fact))
        std::cout << RM << "eigenVectors()] - Hay autovalores duplicados, no se asegura la obtención de todos los autovectores.\n";
    
    for (int i=0; i<n; i++) {
        eVe.Ones();
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
                eVe.Ones();
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
    
    Ainv.Zero();
    eVe.Ones();
    
    eVa = eigenValue+1.0;
    k=1;
    Ainv = (* this) - eigenValue*(1.0+fact)*Identity(n,n);
    Ainv = Ainv^(-1);
    
    while ((eVa < eigenValue*(1.0-fact) || eVa > eigenValue*(1.0+fact)) && k < maxIte) {
        eVe = Ainv * eVe;
        maxS = eVe.maxAbs_sig();
        eVa = 1/maxS + eigenValue*(1.0+fact);
        eVe /= maxS;
        k++;
        
        if ((eVa < eigenValue*(1.0-fact) || eVa > eigenValue*(1.0+fact)) && k == maxIte) {
            eVe.Ones();
            eVe[0][0] *= 2.0;
            eVa = eigenValue+1.0;
            k=1;
            Ainv = (* this) - eigenValue*(1.0+fact)*Identity(n,n);
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
    eVeV = eVe.GetColumnAsVector(0);
    
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
    
    A.Zero();
    Q.Zero();
    R.Zero();
    A = (* this);
    
    while (cota > err && k<maxIte) {
        Q = A.QR(QRmatrix);
        Q||R;
        diagOld = A.GetDiagonal();
        A = R*Q;
        cota = (diagOld - A.GetDiagonal()).AbsoluteMaximumElement();
        err = (A.GetDiagonal()).AbsoluteMaximumElement()*factor;
        k++;
    }
    
    if ((opt& EVa_Ite) != 0)
        std::cout << RM << "eigenValues()] - Iteraciones para calcular los autovalores: " << k << std::endl;
    
    return A.GetDiagonal();
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
template<typename T> inline Matrix<T>
Matrix<T>::operator-()
{
    Matrix<T> tmp(n,m);
    tmp.Zero();
    for (int k=0; k<size; k++) {
        tmp.a[k] = -a[k];
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
    tmp.Zero();
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
Matrix<T>::operator^(const int exp)
{
    Matrix<T> tmp(n,m);
    tmp.Zero();
    if (n != m) {
        std::cout << RM << "operator^exp] - No se puede elevar la matriz a " << exp << " porque no es cuadrada.\n";
        tmp.Ones();
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
operator&&(const Matrix<T>& Mi, const Matrix<T>& Md)
{
    if (Mi.Rows() != Md.Rows())
        std::cout << RM << "operator&&(Mi, Md)] - Estas matrices no tienen la misma altura.\n";
    int n,m;
    n = std::max(Mi.Rows(), Md.Rows());
    m = Mi.Columns() + Md.Columns();
    Matrix<T> tmp(n,m);
    tmp.Zero();
    
    for (int i=0; i<Mi.Rows(); i++) {
        for (int j=0; j<Mi.Columns(); j++) {
            tmp[i][j] = Mi[i][j];
        }
    }
    
    for (int i=0; i<Md.Rows(); i++) {
        for (int j=Mi.Columns(); j<m; j++) {
            tmp[i][j] = Md[i][j-Mi.Columns()];
        }
    }
    
    return tmp;
}

template<typename T> inline void
operator||(Matrix<T>& Mi, Matrix<T>& Md)
{
    if (Mi.Rows() != Md.Rows())
        std::cout << RM << "operator||(Mi, Md)] - Las matrices no tienen la misma altura.\n";
    
    Md = Mi.GetMatrix(0, (Mi.Columns() - Md.Columns()), Md.Rows(), Md.Columns());
    Mi = Mi.GetMatrix(0, 0, Mi.Rows(), (Mi.Columns() - Md.Columns()));
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
    
    M.Resize(n, m);
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
//        if (M.Rows() == 1 && M.Columns() == 1) {
//            out << M(0) << endl;
//        } else if (M.Rows() == 1) {
//            out << "[" << M(0);
//            for (int j=1; j<M.Columns(); j++) {
//                out << "\t" << M(j); 
//            }
//            out << "]\n";
//        } else {
//            out << "⎡" << M(0);
//            for (int j=1; j<M.Columns(); j++) {
//                out << right << "\t" << M(j);
//            }
//            out.width();
//            out << "\t⎤\n";
//            for (int i=1; i<M.Rows()-1; i++) {
//                out << "⎢" << M(i*M.Columns());
//                for (int j=1; j<M.Columns(); j++) {
//                    out << "\t" << M(i*M.Columns() + j);
//                }
//                out << "\t⎥\n";
//            }
//            out << "⎣" << M((M.Rows()-1)*M.Columns());
//            for (int j=1; j<M.Columns(); j++) {
//                out << "\t" << M((M.Rows()-1)*M.Columns() + j);
//            }
//            out << "\t⎦\n";
//        }
//        out.precision();
//    } else {
//        for (int i=0; i<M.Rows(); i++) {
//            for (int j=0; j<M.Columns()-1; j++) {
//                out << M(i*M.Columns() + j) << ",";
//            }
//            out << M(i*M.Columns() + (M.Columns()-1));
//            if (i != M.Rows()-1)
//                out << endl; 
//        }
//    }
    return out;
}



//  --- OTHER FUNCTIONS ---


template<typename T> inline Matrix<T>
cda::math::containers::setdiff(Matrix<T>& A, const Matrix<T>& B, const int& reps)
{
    Matrix<T> tmp(A.Rows(), A.Columns());
    int dim = 0, coincidence = 0;
    
    for (int i=0; i<A.Rows(); i++) {
        
        bool equal = false;
        for (int j=0; j<B.Rows(); j++) {
            
            int count=0;
            for (int k=0; k<A.Columns(); k++) {
                if (A(i,k) == B(j,k))
                    count++;
            }
            
            if (count == A.Columns()) {
                equal = true;
                coincidence++;
                break;
            }
        }
        
        if (equal != true) {
            tmp.setRowV(dim,A.GetRowAsVector(i));
            dim++;
        }
        
        if (coincidence == reps) {
            tmp.setMatrix(dim,0,A.getMatrix(i+1,0));
            dim = A.Rows() - reps;
            break;
        }
    }
    
    if (dim == 0) {
        Matrix<T> fail(1,1);
        fail(0,0) = 0.0;
        return fail;
    } else {
        tmp.Resize(dim, A.Columns());
        return tmp;
    }
}

template<typename T> inline Matrix<T>
cda::math::containers::setdiff(Matrix<T>& A, const Matrix<T>& B)
{
    return setdiff(A, B, A.Rows());
}
