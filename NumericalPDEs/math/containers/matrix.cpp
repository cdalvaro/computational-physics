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

template<typename T> inline Vector<T>
Matrix<T>::solveLU(const Vector<T>& B) const
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
    
    Matrix<T> I(n,n);
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
    
    auto Q(I - (v*vt) * 2.0/(vt.SquaredNorm()));
    auto R(Q * (* this));
    for (int i=1; i<n; i++) {
        R.a[i*m] = (T)0.0;
    }
    
    //  Other columns
    for (int i=1; i<n-1; i++) {
        Matrix<T> H(Q.Rows(), Q.Columns());
        v.Resize(n-i, 1);
        
        c = R.GetMatrix(i, i).GetColumnAsVector(0);
        
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
            auto A(I);
            A.SetMatrix(i, i, H);
            H = A;
        }
        
        Q = H*Q;
        R = H*R;
        
        for (int j=0; j<i; j++) {
            R[i][j] = (T)0.0;
        }
    }
    
    Q = Q.Transpose();
    
    for (int i=1; i<n; i++) {
        for (int j=0; j<i; j++) {
            R[i][j] = (T)0.0;
        }
    }
    
    if (QorR & Qmatrix) {
        return Q;
    } else if (QorR & Rmatrix) {
        return R;
    }
    
    return Q && R;
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
    Ainv = Ainv.Pow(-1);
    
    while ((eVa < eigenValue*(1.0-fact) || eVa > eigenValue*(1.0+fact)) && k < maxIte) {
        eVe = Ainv * eVe;
        maxS = eVe.AbsoluteMaximumElementWithSign();
        eVa = 1/maxS + eigenValue*(1.0+fact);
        eVe /= maxS;
        k++;
        
        if ((eVa < eigenValue*(1.0-fact) || eVa > eigenValue*(1.0+fact)) && k == maxIte) {
            eVe.Ones();
            eVe[0][0] *= 2.0;
            eVa = eigenValue+1.0;
            k=1;
            Ainv = (* this) - eigenValue*(1.0+fact)*Identity(n,n);
            Ainv = Ainv.Pow(-1);
            
            while (1/eVa != eigenValue && k < maxIte) {
                eVe = Ainv * eVe;
                maxS = eVe.AbsoluteMaximumElementWithSign();
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
    
    T cota = 1.0;
    T err = 0.1;
    
    Matrix<T> R(n,n);
    auto A(*this);
    
    size_t k;
    for (k = 1; k <= maxIte; ++k) {
        auto Q(A.QR(QRmatrix));
        Q||R;
        
        auto diagOld(A.GetDiagonal());
        A = R*Q;
        
        cota = (diagOld - A.GetDiagonal()).AbsoluteMaximumElement();
        err = (A.GetDiagonal()).AbsoluteMaximumElement()*factor;
        
        if (cota <= err) {
            break;
        }
    }
    
    if (opt & EVa_Ite) {
        std::cout << RM << "eigenValues()] - Iteraciones para calcular los autovalores: " << k << std::endl;
    }
    
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
