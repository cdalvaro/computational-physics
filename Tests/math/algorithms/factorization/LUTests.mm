//
//  LUTests.m
//  Tests
//
//  Created by Carlos David on 16/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#import <XCTest/XCTest.h>

#import "../../../TestsTools.h"
#import "../../../../computational-physics/math/algorithms/factorization/lu.hpp"
#import "../../../../computational-physics/math/containers/matrix.hpp"

using namespace cda::math::containers;
using namespace cda::math::algorithms::factorization;


@interface LUTests : XCTestCase

@end

@implementation LUTests

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
    [TestsTools setDefaultWorkingDirectory];
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

- (void)testConstructor {
    Matrix<double> matrix_test({
        {3,  2,  1,  2},
        {7,  6,  5,  1},
        {12, 10,  9,  8},
        {15, 14, 13, 12}
    });
    
    XCTAssertNoThrow(LU<Matrix>(matrix_test), "LU matrix constructor OK");
    XCTAssertThrows(LU<Matrix>(Matrix<double>(4, 3)), "LU matrix is not valid for non-square matrices");
}

- (void)testLUMatrices {
    Matrix<int> matrix_test({
        {3,  2,  1,  2},
        {7,  6,  5,  1},
        {12, 10,  9,  8},
        {15, 14, 13, 12}
    });
    
    LU<Matrix, double> lu(matrix_test);
    
    const auto accuracy = 1E-13;
    
    // Matrix U
    Matrix<double> expected_u({
        {3.0,     2.0,     1.0,       2.0},
        {0.0, 4.0/3.0, 8.0/3.0, -11.0/3.0},
        {0.0,     0.0,     1.0,  11.0/2.0},
        {0.0,     0.0,     0.0,      13.0}
    });
    
    XCTAssert([TestsTools compareMatrix:lu.U()
                           withExpected:expected_u
                           whitAccuracy:accuracy],
              "Matrix U OK");
    
    // Matrix L
    Matrix<double> expected_l({
        {    1.0,     0.0, 0.0, 0.0},
        {7.0/3.0,     1.0, 0.0, 0.0},
        {    4.0, 3.0/2.0, 1.0, 0.0},
        {    5.0,     3.0, 0.0, 1.0}
    });
    
    XCTAssert([TestsTools compareMatrix:lu.L()
                           withExpected:expected_l
                           whitAccuracy:accuracy],
              "Matrix L OK");
}

- (void)testDeterminant {
    const Matrix<double> matrix1({
        { 21, 18, 15,  4},
        { 49, 41, 35,  7},
        { 84, 72, 63, 12},
        {105, 90, 75, 15}
    });
    
    LU<Matrix, double> lu1(matrix1);
    XCTAssertEqual(lu1.determinant(), 315, "determinant OK");
    
    const Matrix<double> matrix2({
        { 21, 18, 15,  4},
        { 42, 36, 30,  8},
        { 84, 72, 63, 12},
        {105, 90, 75, 15}
    });
    
    LU<Matrix, double> lu2(matrix2);
    XCTAssertEqual(lu2.determinant(), 0, "determinant 0 OK");
}

- (void)testInverseMatrix {
    const Matrix<double> matrix1({
        { 3,  2,  4},
        { 7,  6,  5},
        {11, 10,  9}
    });
    
    const auto expected = Matrix<double>({
        { 4,  22, -14},
        {-8, -17,  13},
        { 4,  -8,   4}
    }) / 12.0;
    
    LU<Matrix, double> lu1(matrix1);
    XCTAssert([TestsTools compareMatrix:lu1.InverseMatrix()
                           withExpected:expected
                           whitAccuracy:TESTS_TOOLS_DEFAULT_ACCURACY],
              "InverseMatrix OK");
    
    const Matrix<double> matrix2({
        {0,  1,  2,  3},
        { 4,  5,  6,  7},
        { 8,  10, 7, 14},
        {12, 13, 14, 15}
    });
    
    LU<Matrix, double> lu2(matrix2);
    XCTAssertThrows(lu2.InverseMatrix(), "Matrix is degenerate");
}

- (void)testSolveLinearSystem {
    const Matrix<double> matrix({
        { 1,  0, 1},
        { 0, -3, 1},
        { 2,  1, 3}
    });
    
    LU<Matrix, double>lu(matrix);
    
    const Vector<double> terms({6, 7, 15});
    
    const Vector<double> expected({2, -1, 4});
    XCTAssert([TestsTools compareVector:lu.SolveLinearSystem(terms)
                           withExpected:expected
                           whitAccuracy:TESTS_TOOLS_DEFAULT_ACCURACY],
              "SolveLinearSystem OK");
    
    XCTAssertThrows(lu.SolveLinearSystem(Vector<double>({1, 2, 3, 4})),
                    "The number of independent terms does not match the matrix dimension");
}

@end
