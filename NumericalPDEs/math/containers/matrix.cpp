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






//  --- MÉTODO QR ---

//template<typename T> inline Matrix<T>
//Matrix<T>::eigenVectors(int maxIte, unsigned char opt)
//{
//    if (n != m)
//        std::cout << RM << "eigenVectors()] -> ";
//    
//    int k;
//    T eVa, maxS = 0.0, fact = 1E-07;
//    Vector<T> eVas(n);
//    Matrix<T> eVes(n,n), eVe(n,1), Ainv(n,n);
//    
//    Ainv.Zero();
//    eVas.Zero();
//    eVes.Zero();
//    
//    eVas = this -> eigenValues(maxIte, preDef, opt);
//    if ((opt& ORGANIZED) != 0)
//        eVas = eVas.organize();
//    
//    if (eVas.duplicate(fact))
//        std::cout << RM << "eigenVectors()] - Hay autovalores duplicados, no se asegura la obtención de todos los autovectores.\n";
//    
//    for (int i=0; i<n; i++) {
//        eVe.Ones();
//        eVa = (T)eVas(i)+1.0;
//        k=1;
//        Ainv = (* this) - identity(n,n)*eVas(i)*(1.0+fact);
//        Ainv = Ainv^(-1);
//        
//        while ((eVa < eVas(i)*(1.0-fact) || eVa > eVas(i)*(1.0+fact)) && k < maxIte) {
//            eVe = Ainv * eVe;
//            maxS = eVe.maxAbs_sig();
//            eVa = 1/maxS + eVas(i)*(1.0+fact);
//            eVe /= maxS;
//            k++;
//            
//            if ((eVa < eVas(i)*(1.0-fact) || eVa > eVas(i)*(1.0+fact)) && k == maxIte) {
//                eVe.Ones();
//                eVe(i) *= 2.0;
//                eVa = (T)eVas(i)+1.0;
//                k=1;
//                Ainv = (* this) - identity(n,n)*eVas(i)*(1.0+fact);
//                Ainv = Ainv^(-1);
//                
//                while (1/eVa != eVas(i) && k < maxIte) {
//                    eVe = Ainv * eVe;
//                    maxS = eVe.maxAbs_sig();
//                    eVa = 1/maxS + eVas(i)*(1.0+fact);
//                    eVe /= maxS;
//                    k++;
//                }
//            }
//        }
//        eVes.setColumn(i, eVe);
//        
//        if ((opt& EVe_Ite) != 0)
//            std::cout << RM << "eigenVectors()] - Iteraciones para calcular el autovector " << i+1 << ": " << k << std::endl;
//    }
//    return eVes;
//}


