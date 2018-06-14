//
//  MatrixTests.m
//  Tests
//
//  Created by Carlos David on 11/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#import <XCTest/XCTest.h>

#import "../../../NumericalPDEs/math/containers/matrix.hpp"

using namespace cda::math::containers;

@interface MatrixTests : XCTestCase

@end

@implementation MatrixTests

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

- (void)testMatrixProduct {
    const Matrix<int> matrix1(4, 4, {
         0,  1,  2,  3,
         4,  5,  6,  7,
         8,  9, 10, 11,
        12, 13, 14, 15
    });
    
    const Matrix<int> matrix2(4, 4, {
         3,  2,  1,  0,
         7,  6,  5,  4,
        11, 10,  9,  8,
        15, 14, 13, 12
    });
    
    const Matrix<int> expected(4, 4, {
        74, 68, 62, 56,
        218, 196, 174, 152,
        362, 324, 286, 248,
        506, 452, 398, 344
    });
    
    const auto result = matrix1 * matrix2;
    
    XCTAssert(result == expected, "Matrix product OK");
}

- (void)testTransposeOperation {
    const Matrix<int> matrix({
        {0,  1,  2,  3,  4},
        {5,  6,  7,  8,  9},
        {10, 11, 12, 13, 14},
        {15, 16, 17, 18, 19}
    });
    
    const Matrix<int> expected({
        {0, 5, 10, 15},
        {1, 6, 11, 16},
        {2, 7, 12, 17},
        {3, 8, 13, 18},
        {4, 9, 14, 19}
    });
    
    const auto result = matrix.Transpose();
    
    XCTAssert(result == expected, "Transpose operation OK");
}

- (void)testGetMatrixMethods {
    const Matrix<int> matrix({
        {0,  1,  2,  3,  4},
        {5,  6,  7,  8,  9},
        {10, 11, 12, 13, 14},
        {15, 16, 17, 18, 19}
    });
    
    const Matrix<int> expected1({
        {7,  8,  9},
        {12, 13, 14},
        {17, 18, 19}
    });
    
    const auto result1 = matrix.GetMatrix(1, 2);
    
    XCTAssert(result1 == expected1, "GetMatrix without lengths OK");
    
    const Matrix<int> expected2({
        {0,  1,  2},
        {5,  6,  7},
        {10, 11, 12}
    });
    
    const auto result2 = matrix.GetMatrix(0, 0, 3, 3);
    
    XCTAssert(result2 == expected2, "GetMatrix with lengths OK");
}

- (void)testSetMatrixMethods {
    Matrix<int> result(4, 5, 0);
    result.SetMatrix(1, 2,
                     Matrix<int>({
        {7,  8,  9},
        {12, 13, 14},
        {17, 18, 19}
    }));
    
    const Matrix<int> expected1({
        {0, 0,  0,  0,  0},
        {0, 0,  7,  8,  9},
        {0, 0, 12, 13, 14},
        {0, 0, 17, 18, 19}
    });
    
    XCTAssert(result == expected1, "SetMatrix without lengths OK");
    
    result.SetMatrix(0, 0,
                     Matrix<int>({
        { 0,  1,  2},
        { 5,  6,  7},
        {10, 11, 12}
    }));
    
    const Matrix<int> expected2({
        { 0,  1,  2,  0,  0},
        { 5,  6,  7,  8,  9},
        {10, 11, 12, 13, 14},
        { 0,  0, 17, 18, 19}
    });
    
    XCTAssert(result == expected2, "SetMatrix with lengths OK");
}

- (void)testMarrixLUMethods {
    Matrix<double> matrix_test({
        {3,  2,  1,  2},
        {7,  6,  5,  1},
        {12, 10,  9,  8},
        {15, 14, 13, 12}
    });
    
    const auto accuracy = 1E-13;
    
    // Matrix U
    Matrix<double> expected_u({
        {3.0,     2.0,     1.0,       2.0},
        {0.0, 4.0/3.0, 8.0/3.0, -11.0/3.0},
        {0.0,     0.0,     1.0,  11.0/2.0},
        {0.0,     0.0,     0.0,      13.0}
    });
    
    auto result_u = matrix_test.U();
    
    for (auto it_result = result_u.Begin(), it_expected = expected_u.Begin();
         it_result != result_u.End(); ++it_result, ++it_expected) {
        XCTAssertEqualWithAccuracy(*it_result, *it_expected, accuracy, "Element of matrix U is equal");
    }
    
    // Matrix L
    Matrix<double> expected_l({
        {    1.0,     0.0, 0.0, 0.0},
        {7.0/3.0,     1.0, 0.0, 0.0},
        {    4.0, 3.0/2.0, 1.0, 0.0},
        {    5.0,     3.0, 0.0, 1.0}
    });
    
    auto result_l = matrix_test.L();
    
    for (auto it_result = result_l.Begin(), it_expected = expected_l.Begin();
         it_result != result_l.End(); ++it_result, ++it_expected) {
        XCTAssertEqualWithAccuracy(*it_result, *it_expected, accuracy, "Element of matrix L is equal");
    }
    
}

- (void)testSumRows {
    Matrix<int> matrix_test({
        { 3,  2,  1,  2},
        { 7,  6,  5,  1},
        {12, 10,  9,  8},
        {15, 14, 13, 12}
    });
    
    Matrix<int> expected1(4, 1, {8, 19, 39, 54});
    auto result1 = matrix_test.SumRows();
    
    XCTAssert(result1 == expected1, "SumRows OK");
    
    auto expected2 = expected1.GetColumnAsVector(0);
    auto result2 = matrix_test.SumRowsAsVector();
    
    XCTAssert(result2 == expected2, "SumRowsAsVector OK");
}

- (void)testSumColumns {
    Matrix<int> matrix_test({
        { 3,  2,  1,  2},
        { 7,  6,  5,  1},
        {12, 10,  9,  8},
        {15, 14, 13, 12}
    });
    
    Matrix<int> expected1(1, 4, {37, 32, 28, 23});
    auto result1 = matrix_test.SumColumns();
    
    XCTAssert(result1 == expected1, "CumColumns OK");
    
    auto expected2 = expected1.GetRowAsVector(0);
    auto result2 = matrix_test.SumColumnsAsVector();
    
    XCTAssert(result2 == expected2, "SumColumnsAsVector OK");
}

@end
