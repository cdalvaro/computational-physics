//
//  MatrixTests.m
//  Tests
//
//  Created by Carlos David on 11/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#import <XCTest/XCTest.h>

#import "../../TestsTools.h"
#import "../../../computational-physics/math/containers/matrix.hpp"

using namespace cda::math::containers;


@interface MatrixTests : XCTestCase

@end

@implementation MatrixTests

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
    [TestsTools setDefaultWorkingDirectory];
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

- (void)testConstructors {
    
    const Matrix<double> matrix;
    XCTAssert(matrix.is_empty(), "matrix is empty");
    XCTAssert(matrix.is_null(), "matrix is null");
    
    const Matrix<double> matrix1(4, 4, 1);
    XCTAssert(!matrix1.is_empty() && !matrix1.is_null(), "matrix1 is not empty and is not null");
    XCTAssertEqual(matrix1.Rows(), 4, "The number of rows of matrix1 is OK");
    XCTAssertEqual(matrix1.Columns(), 4, "The number of columns of matrix1 is OK");
    XCTAssertEqual(matrix1.size(), 16, "The size of matrix1 is OK");
    
    bool all_elements_are_one = true;
    for (auto it = matrix1.begin(); it != matrix1.end(); ++it) {
        if (*it != 1) {
            all_elements_are_one = false;
            break;
        }
    }
    XCTAssert(all_elements_are_one, "All elements of matrix1 are OK");
    
    const Matrix<double> matrix2(4, 4, {
        0,  1,  2,  3,
        4,  5,  6,  7,
        8,  9, 10, 11,
        12, 13, 14, 15
    });
    XCTAssert(!matrix2.is_empty() && !matrix2.is_null(), "matrix2 is not empty and is not null");
    XCTAssertEqual(matrix2.Rows(), 4, "The number of rows of matrix2 is OK");
    XCTAssertEqual(matrix2.Columns(), 4, "The number of columns of matrix2 is OK");
    XCTAssertEqual(matrix2.size(), 16, "The size of matrix2 is OK");
    
    XCTAssertThrows(Matrix<double>(2, 2, { 0,  1,  2 }), "The array of values has different size than the matrix");
    
    const Matrix<double> matrix3({
        { 0,  1,  2,  3},
        { 4,  5,  6,  7},
        { 8,  9, 10, 11},
        {12, 13, 14, 15}
    });
    XCTAssert(!matrix3.is_empty() && !matrix3.is_null(), "matrix3 is not empty and is not null");
    XCTAssertEqual(matrix3.Rows(), 4, "The number of rows of matrix3 is OK");
    XCTAssertEqual(matrix3.Columns(), 4, "The number of columns of matrix3 is OK");
    XCTAssertEqual(matrix3.size(), 16, "The size of matrix3 is OK");
    
    Matrix<double> matrix4(matrix3);
    XCTAssertEqual(matrix4, matrix3, "matrix3 and matrix4 are equal");
    
    Matrix<double> matrix5(std::move(matrix4));
    XCTAssertEqual(matrix5, matrix3, "matrix3 and matrix4 are equal");
    XCTAssert(matrix4.is_empty(), "matrix4 is empty after had been moved");
    XCTAssert(matrix4.is_null(), "matrix4 is null after had been moved");
}

- (void)testAssign {
    const Matrix<double> matrix(4, 4, {
        0,  1,  2,  3,
        4,  5,  6,  7,
        8,  9, 10, 11,
        12, 13, 14, 15
    });
    
    Matrix<double> test;
    test = matrix;
    
    XCTAssertEqual(test, matrix, "Assignment OK");
    
    test = test;
    XCTAssertEqual(test, matrix, "Self assignment OK");
    
    Matrix<double> test2;
    test2 = std::move(test);
    
    XCTAssertEqual(test2, matrix, "Move assignment OK");
    XCTAssert(test.is_null(), "test matrix is null after had been moved");

//Save the diagnostic state
#pragma clang diagnostic push
    
// Ignore: Explicitly moving variable of type 'Matrix<double>' to itself
#pragma clang diagnostic ignored "-Wself-move"
    test2 = std::move(test2);
    XCTAssertEqual(test2, matrix, "Move assignment does not destroy object if self matrix OK");
    
//Restore the disgnostic state
#pragma clang diagnostic pop
}

- (void)testOnes {
    const Matrix<double> expected_ones(5, 5, 1);
    XCTAssertEqual(Matrix<double>::ones(5, 5), expected_ones, "ones static method OK");
    
    Matrix<double> result(5, 5, 2);
    result.ones();
    XCTAssertEqual(result, expected_ones, "ones method OK");
}

- (void)testZero {
    const Matrix<double> expected_zero(5, 5, 0);
    XCTAssertEqual(Matrix<double>::zero(5, 5), expected_zero, "zero static method OK");
    
    Matrix<double> result(5, 5, 2);
    result.zero();
    XCTAssertEqual(result, expected_zero, "zero method OK");
}

- (void)testIdentity {
    const Matrix<double> expected_identity({
        {1, 0, 0, 0, 0},
        {0, 1, 0, 0, 0},
        {0, 0, 1, 0, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 0, 0, 1}
    });
    XCTAssertEqual(Matrix<double>::Identity(5), expected_identity, "Identity static method OK");
    
    Matrix<double> result(5, 5, 2);
    result.Identity();
    XCTAssertEqual(result, expected_identity, "Identity method OK");
    
    Matrix<double> result2(5, 4, 2);
    XCTAssertThrows(result2.Identity(), "Identity method is only valid for square matrices");
}

- (void)testComparison {
    const Matrix<double> matrix1({
        { 0,  1,  2,  3},
        { 4,  5,  6,  7},
        { 8,  9, 10, 11},
        {12, 13, 14, 15}
    });
    
    XCTAssert(matrix1 == matrix1, "Matrix is equal to itself");
    XCTAssert(!(matrix1 != matrix1), "Matrix is not different to itself");
    
    const Matrix<double> matrix2({
        { 0,  1,  2,  3},
        { 8,  9, 10, 11},
        { 4,  5,  6,  7},
        {12, 13, 14, 15}
    });
    
    XCTAssert(matrix1 != matrix2, "matrix1 is different to matrix2");
    XCTAssert(!(matrix1 == matrix2), "matrix1 is not equal to matrix2");
    
    // Same as matrix1 whit less columns
    const Matrix<double> matrix3({
        { 0,  1,  2},
        { 4,  5,  6},
        { 8,  9, 10},
        {12, 13, 14}
    });
    
    XCTAssert(matrix3 != matrix1, "matrix3 is different to matrix1");
    XCTAssert(!(matrix3 == matrix1), "matrix3 is not equal to matrix1");
    
    // Same as matrix2 whit less rows
    const Matrix<double> matrix4({
        { 0,  1,  2,  3},
        { 4,  5,  6,  7},
        {12, 13, 14, 15}
    });
    
    XCTAssert(matrix4 != matrix2, "matrix4 is different to matrix2");
    XCTAssert(!(matrix4 == matrix2), "matrix4 is not equal to matrix2");
}

- (void)testResize {
    const Matrix<double> matrix({
        {0, 1,  2},
        {4, 5,  6},
        {8, 9, 10}
    });
    
    auto result(matrix);
    result.resize(matrix.Rows(), matrix.Columns());
    XCTAssertEqual(result, matrix, "matrix has not been changed");
    
    const Matrix<double> expected1({
        {0, 1},
        {4, 5}
    });
    
    result = Matrix<double>(matrix);
    result.resize(2, 2);
    XCTAssertEqual(result, expected1, "matrix has been resized OK to a 2x2 matrix");
    
    const Matrix<double> expected2({
        {0, 1},
        {4, 5},
        {8, 9}
    });
    
    result = Matrix<double>(matrix);
    result.resize(3, 2);
    XCTAssertEqual(result, expected2, "matrix has been resized OK to a 3x2 matrix");
    
    const Matrix<double> expected3({
        {0, 1, 2},
        {4, 5, 6}
    });
    
    result = Matrix<double>(matrix);
    result.resize(2, 3);
    XCTAssertEqual(result, expected3, "matrix has been resized OK to a 2x3 matrix");
    
    const Matrix<double> expected4({
        {0, 1,  2, 0},
        {4, 5,  6, 0},
        {8, 9, 10, 0},
        {0, 0,  0, 0}
    });
    
    result = Matrix<double>(matrix);
    result.resize(4, 4, true);
    XCTAssertEqual(result, expected4, "matrix has been resized OK to a 4x4 matrix adding zeros to the new elements");
    
    const Matrix<double> expected5({
        {0, 1,  2},
        {4, 5,  6},
        {8, 9, 10},
        {0, 0,  0}
    });
    
    result = Matrix<double>(matrix);
    result.resize(4, 3, true);
    XCTAssertEqual(result, expected5, "matrix has been resized OK to a 4x3 matrix adding zeros to the new elements");
}

- (void)testChangeDimensions {
    const Matrix<double> matrix({
        {0, 1,  2, 4},
        {4, 5,  6, 7},
        {8, 9, 10, 11}
    });
    
    const Matrix<double> expected({
        {0,  1,  2},
        {4,  4,  5},
        {6,  7,  8},
        {9, 10, 11}
    });
    
    auto result(matrix);
    result.change_dimensions(4, 3);
    
    XCTAssertEqual(result, expected, "Dimensions have been changed OK");
    XCTAssertThrows(result.change_dimensions(3, 3), "The total size of the matrix cannot be changed");
}

- (void)testProducts {
    const Matrix<double> matrix1({
        { 0,  1,  2,  3},
        { 4,  5,  6,  7},
        { 8,  9, 10, 11},
        {12, 13, 14, 15}
    });
    
    const Matrix<double> matrix2({
        { 3,  2,  1},
        { 7,  6,  5},
        {11, 10,  9},
        {15, 14, 13}
    });
    
    const Matrix<double> expected1({
        { 74,  68,  62},
        {218, 196, 174},
        {362, 324, 286},
        {506, 452, 398}
    });
    
    XCTAssertEqual(matrix1 * matrix2, expected1, "Square matrices product OK");
    XCTAssertThrows(matrix2 * matrix1, "Incompatible dimensions to compute the product");
    
    const Matrix<double> matrix3({
        { 0,  1,  2,  3, -1},
        { 4,  5,  6,  7, -2},
        { 8,  9, 10, 11, -3},
        {12, 13, 14, 15, -4}
    });
    
    const Matrix<double> matrix4({
        { 3,   2,   1,   0, -11, 235},
        { 7,   6,   5,   4, -22, 264},
        {11,  10,   9,   8, -33, 436},
        {15,  14,  13,  12, -44, 643},
        {54, 235,  21,  21, 235, 426}
    });
    
    const Matrix<double> expected2({
        { 20, -167,  41,  35,  -455,  2639},
        {110, -274, 132, 110, -1130,  8525},
        {200, -381, 223, 185, -1805, 14411},
        {290, -488, 314, 260, -2480, 20297}
    });
    
    XCTAssertEqual(matrix3 * matrix4, expected2, "Rectangular matrices product OK");
    
    XCTAssertThrows(matrix4 * matrix3, "Matrices dimensions are incompatible");
    
    const Matrix<double> expected3({
        { 0,  2,  4,  6},
        { 8, 10, 12, 14},
        {16, 18, 20, 22},
        {24, 26, 28, 30}
    });
    
    XCTAssertEqual(matrix1 * 2, expected3, "Product between matrix and scalar OK");
    XCTAssertEqual(2.0 * matrix1, expected3, "Product between scalar and matrix OK");
    
    Matrix<double> matrix1_copy({
        { 0,  1,  2,  3},
        { 4,  5,  6,  7},
        { 8,  9, 10, 11},
        {12, 13, 14, 15}
    });
    
    XCTAssertEqual(matrix1_copy *= 2, expected3, "Product between matrix and scalar OK");
}

- (void)testDivisions {
    Matrix<double> matrix1({
        { 0,  1,  2,  3},
        { 4,  5,  6,  7},
        { 8,  9, 10, 11},
        {12, 13, 14, 15}
    });
    
    Matrix<double> expected({
        {0, 0.5, 1, 1.5},
        {2, 2.5, 3, 3.5},
        {4, 4.5, 5, 5.5},
        {6, 6.5, 7, 7.5}
    });
    
    XCTAssertEqual(matrix1 / 2.0, expected, "Division OK");
    XCTAssertEqual(matrix1 /= 2.0, expected, "Division OK");
}

- (void)testAdditionOfMatrices {
    Matrix<double> matrix1({
        { 3,  2,  1},
        { 7,  6,  5},
        {11, 10,  9}
    });
    
    const Matrix<double> matrix2({
        { 4,  62,    6},
        {34,  73,  375},
        {25, 251, 2531}
    });
    
    const Matrix<double> expected({
        { 7,  64,    7},
        {41,  79,  380},
        {36, 261, 2540}
    });
    
    XCTAssertEqual(matrix1 + matrix2, expected, "Addition of matrix1 plus matrix2 OK");
    XCTAssertEqual(matrix2 + matrix1, expected, "Addition of matrix2 plus matrix1 OK");
    XCTAssertEqual(matrix1 += matrix2, expected, "Addition of matrix2 over matrix1 OK");
    
    const Matrix<double> matrix3({
        { 4,  62,    6},
        {34,  73,  375},
        {25, 251, 2531},
        {34, 215,  321}
    });
    
    const Matrix<double> matrix4({
        { 4,  62,    6,  132},
        {34,  73,  375, 3215},
        {25, 251, 2531, 3125}
    });
    
    XCTAssertThrows(matrix2 + matrix3, "Add incompatible dimensions");
    XCTAssertThrows(matrix2 + matrix4, "Add incompatible dimensions");
    XCTAssertThrows(matrix1 += matrix3, "Add incompatible dimensions");
    XCTAssertThrows(matrix1 += matrix3, "Add incompatible dimensions");
}

- (void)testNegativeMatrix {
    const Matrix<double> matrix(4, 4, 1);
    XCTAssertEqual(-matrix, matrix * -1, "Negative matrix OK");
}

- (void)testSubtractionOfMatrices {
    Matrix<double> matrix1({
        { 3,  2,  1},
        { 7,  6,  5},
        {11, 10,  9}
    });
    
    const Matrix<double> matrix2({
        { 4,  62,    6},
        {34,  73,  375},
        {25, 251, 2531}
    });
    
    const Matrix<double> expected({
        { -1,  -60,    -5},
        {-27,  -67,  -370},
        {-14, -241, -2522}
    });
    
    XCTAssertEqual(matrix1 - matrix2, expected, "Subtraction of matrix1 plus matrix2 OK");
    XCTAssertEqual(matrix2 - matrix1, -expected, "Subtraction of matrix2 plus matrix1 OK");
    XCTAssertEqual(matrix1 -= matrix2, expected, "Subtraction of matrix2 over matrix1 OK");
    
    const Matrix<double> matrix3({
        { 4,  62,    6},
        {34,  73,  375},
        {25, 251, 2531},
        {34, 215,  321}
    });
    
    const Matrix<double> matrix4({
        { 4,  62,    6,  132},
        {34,  73,  375, 3215},
        {25, 251, 2531, 3125}
    });
    
    XCTAssertThrows(matrix2 - matrix3, "Subtract incompatible dimensions");
    XCTAssertThrows(matrix2 - matrix4, "Subtract incompatible dimensions");
    XCTAssertThrows(matrix1 -= matrix3, "Subtract incompatible dimensions");
    XCTAssertThrows(matrix1 -= matrix3, "Subtract incompatible dimensions");
}

- (void)testPowers {
    const Matrix<double> matrix1({
        {-1,  1,  2,  3},
        { 4,  5,  6,  7},
        { 8,  10, 7, 14},
        {12, 13, 14, 15}
    });
    
    const auto expected1 = Matrix<double>::Identity(4);
    
    XCTAssertEqual(matrix1.Pow(0), expected1, "Power 0 OK");
    
    const Matrix<double> expected2({
        { 57,  63,  60,  77},
        {148, 180, 178, 236},
        {256, 310, 321, 402},
        {332, 412, 410, 548}
    });
    
    XCTAssertEqual(matrix1.Pow(2), expected2, "Power 2 OK");
    
    const Matrix<double> expected3({
        { 1599,  1973,  1990,  2607},
        { 4828,  5896,  5926,  7736},
        { 8376, 10242, 10247, 13462},
        {11172, 13616, 13678, 17840}
    });
    
    XCTAssertEqual(matrix1.Pow(3), expected3, "Power 3 OK");
    
    const Matrix<double> matrix2({
        { 3,  2,  4},
        { 7,  6,  5},
        {11, 10,  9}
    });
    
    const auto expected4 = Matrix<double>({
        { 4,  22, -14},
        {-8, -17,  13},
        { 4,  -8,   4}
    }) / 12.0;
    
    XCTAssert([TestsTools compareMatrix:matrix2.Pow(-1)
                           withExpected:expected4
                           whitAccuracy:TESTS_TOOLS_DEFAULT_ACCURACY],
              "Power -1 OK");
    
    const Matrix<double> matrix3({
        {0,  1,  2,  3},
        { 4,  5,  6,  7},
        { 8,  10, 7, 14},
        {12, 13, 14, 15}
    });
    XCTAssertThrows(matrix3.Pow(-1), "Matrix is degenerate");
    
    const Matrix<double> matrix4({
        { -1,  1,  2,  3},
        { 4,  5,  6,  7},
        { 8,  9, 10, 11}
    });
    
    XCTAssertThrows(matrix4.Pow(0), "Matrix must be square");
    
    const Matrix<double> matrix5({
        {-1,  1,  2,  3},
        { 4,  5,  6,  7},
        { 8, 10, 12, 14},   // Same row as above multiplied by 2
        {12, 13, 14, 15}
    });
    
    XCTAssertThrows(matrix5.Pow(0), "Matrix is singular");
}

- (void)testTransposeOperation {
    const Matrix<double> matrix({
        {0,  1,  2,  3,  4},
        {5,  6,  7,  8,  9},
        {10, 11, 12, 13, 14},
        {15, 16, 17, 18, 19}
    });
    
    const Matrix<double> expected({
        {0, 5, 10, 15},
        {1, 6, 11, 16},
        {2, 7, 12, 17},
        {3, 8, 13, 18},
        {4, 9, 14, 19}
    });
    
    XCTAssertEqual(matrix.Transpose(), expected, "Transpose operation OK");
}

- (void)testSquareMatrix {
    const auto square_matrix = Matrix<double>::ones(4, 4);
    XCTAssert(square_matrix.IsSquare(), "square_matrix is square");
    
    const auto rectangular_matrix = Matrix<double>::ones(4, 3);
    XCTAssert(!rectangular_matrix.IsSquare(), "rectangular_matrix is not square");
}

- (void)testGettersMethods {
    const Matrix<double> matrix({
        {0,  1,  2,  3,  4},
        {5,  6,  7,  8,  9},
        {10, 11, 12, 13, 14},
        {15, 16, 17, 18, 19}
    });
    
    // --- get_matrix
    const Matrix<double> expected_matrix1({
        {7,  8,  9},
        {12, 13, 14},
        {17, 18, 19}
    });
    
    XCTAssertEqual(matrix.get_matrix(1, 2), expected_matrix1, "get_matrix without lengths OK");
    
    const Matrix<double> expected_matrix2({
        {0,  1,  2},
        {5,  6,  7},
        {10, 11, 12}
    });
    
    XCTAssertEqual(matrix.get_matrix(0, 0, 3, 3), expected_matrix2, "get_matrix with lengths OK");
    
    XCTAssertThrows(matrix.get_matrix(1, 1, 6, 6), "Unable to get submatrix. Elements out of bounds");
    
    // --- get_row
    const Matrix<double> expected_first_row({{0,  1,  2,  3,  4}});
    XCTAssertEqual(matrix.get_row(0), expected_first_row, "get_row first OK");
    
    const Matrix<double> expected_last_row({{15, 16, 17, 18, 19}});
    XCTAssertEqual(matrix.get_row(3), expected_last_row, "get_row last OK");
    
    const Matrix<double> expected_in_between_row({{10, 11, 12, 13, 14}});
    XCTAssertEqual(matrix.get_row(2), expected_in_between_row, "get_row last OK");
    
    XCTAssertThrows(matrix.get_row(4), "get_row out of bounds");
    
    // --- get_column
    const Matrix<double> expected_first_column(4, 1, {0, 5, 10, 15});
    XCTAssertEqual(matrix.get_column(0), expected_first_column, "get_row first OK");
    
    const Matrix<double> expected_last_column(4, 1, {4, 9, 14, 19});
    XCTAssertEqual(matrix.get_column(4), expected_last_column, "get_row last OK");
    
    const Matrix<double> expected_in_between_column(4, 1, {2, 7, 12, 17});
    XCTAssertEqual(matrix.get_column(2), expected_in_between_column, "get_row last OK");
    
    XCTAssertThrows(matrix.get_column(5), "get_column out of bounds");
    
    // --- get_diagonal
    XCTAssertThrows(matrix.get_diagonal(), "This method is only available for square matrices");
    
    const Matrix<double> square_matrix({
        { 0,  1,  2,  3},
        { 5,  6,  7,  8},
        {10, 11, 12, 13},
        {15, 16, 17, 18}
    });
    
    const Vector<double> expected_diagonal({0, 6, 12, 18});
    XCTAssertEqual(square_matrix.get_diagonal(), expected_diagonal, "get_diagonal method OK");
}

- (void)testGettersAsVectorMethods {
    const Matrix<double> matrix({
        {0,  1,  2,  3,  4},
        {5,  6,  7,  8,  9},
        {10, 11, 12, 13, 14},
        {15, 16, 17, 18, 19}
    });
    
    // --- get_column_as_vector ---
    // Get the whole first column
    Vector<double> expected_vector({0, 5, 10, 15});
    XCTAssertEqual(matrix.get_column_as_vector(0), expected_vector, "get_column_as_vector for the whole first column OK!");
    
    // Get the last two elements of the first column
    expected_vector = Vector<double>({10, 15});
    XCTAssertEqual(matrix.get_column_as_vector(0, 2), expected_vector, "get_column_as_vector for the last two elements of the first column OK!");
    
    // Get the whole third column
    expected_vector = Vector<double>({2, 7, 12, 17});
    XCTAssertEqual(matrix.get_column_as_vector(2), expected_vector, "get_column_as_vector for the whole third column OK!");
    
    // Get the last three elements of the fourth column
    expected_vector = Vector<double>({8, 13, 18});
    XCTAssertEqual(matrix.get_column_as_vector(3, 1), expected_vector, "get_column_as_vector for the last three elements of the fourth column OK!");
    
    // Get the last column
    expected_vector = Vector<double>({4, 9, 14, 19});
    XCTAssertEqual(matrix.get_column_as_vector(4), expected_vector, "get_column_as_vector for the last column OK!");
    
    // Get column out of bounds
    XCTAssertThrows(matrix.get_column_as_vector(3, 6), "get_column_as_vector out of bounds by number of elements");
    XCTAssertThrows(matrix.get_column_as_vector(6), "get_column_as_vector out of bounds by number of column");
    
    // --- get_row_as_vector ---
    // Get the whole first row
    expected_vector = Vector<double>({0,  1,  2,  3,  4});
    XCTAssertEqual(matrix.get_row_as_vector(0), expected_vector, "get_row_as_vector for the whole first row OK!");
    
    // Get the last two elements of the first column
    expected_vector = Vector<double>({7,  8,  9});
    XCTAssertEqual(matrix.get_row_as_vector(1, 2), expected_vector, "get_row_as_vector for the last three elements of the second row OK!");
    
    // Get the whole third row
    expected_vector = Vector<double>({10, 11, 12, 13, 14});
    XCTAssertEqual(matrix.get_row_as_vector(2), expected_vector, "get_row_as_vector for the whole third row OK!");
    
    // Get the last element of the last row
    expected_vector = Vector<double>({19});
    XCTAssertEqual(matrix.get_row_as_vector(3, 4), expected_vector, "get_row_as_vector for the last element of the last row OK!");
    
    // Get row out of bounds
    XCTAssertThrows(matrix.get_row_as_vector(3, 6), "get_column_as_vector out of bounds by number of elements");
    XCTAssertThrows(matrix.get_row_as_vector(6), "get_column_as_vector out of bounds by number of column");
}

- (void)testSetMatrixMethods {
    Matrix<double> result(4, 5, 0);
    result.SetMatrix(1, 2,
                     Matrix<double>({
        {7,  8,  9},
        {12, 13, 14},
        {17, 18, 19}
    }));
    
    const Matrix<double> expected1({
        {0, 0,  0,  0,  0},
        {0, 0,  7,  8,  9},
        {0, 0, 12, 13, 14},
        {0, 0, 17, 18, 19}
    });
    
    XCTAssertEqual(result, expected1, "SetMatrix without lengths OK");
    
    result.SetMatrix(0, 0,
                     Matrix<double>({
        { 0,  1,  2},
        { 5,  6,  7},
        {10, 11, 12}
    }));
    
    const Matrix<double> expected2({
        { 0,  1,  2,  0,  0},
        { 5,  6,  7,  8,  9},
        {10, 11, 12, 13, 14},
        { 0,  0, 17, 18, 19}
    });
    
    XCTAssertEqual(result, expected2, "SetMatrix with lengths OK");
    
    result.SetMatrix(0, 0,
                     Matrix<double>({
        { 0,  1,  2,  0,  0, 346},
        { 5,  6,  7,  8,  9,  21},
        {10, 11, 12, 13, 14,   1},
        { 0,  0, 17, 18, 19, 251},
        { 5,  6,  7,  8,  9,  21}
    }));
    
    XCTAssertEqual(result, expected2, "SetMatrix with greater lengths OK");
    
    XCTAssertThrows(result.SetMatrix(5, 4, expected2), "SetMatrix out of bounds by row");
    XCTAssertThrows(result.SetMatrix(3, 6, expected2), "SetMatrix out of bounds by column");
}

- (void)testSetRow {
    auto result = Matrix<double>::ones(4, 5);
    
    const Matrix<double> expected({
        { 1,  1, 1, 1,  1},
        { 2, -5, 7, 9, -7},
        { 1,  1, 1, 1,  1},
        { 1,  1, 1, 1,  1}
    });
    
    Matrix<double> new_row({{2, -5, 7, 9, -7}});
    result.set_row(1, new_row);
    XCTAssertEqual(result, expected, "set_row complete OK");
    
    XCTAssertThrows(result.set_row(4, new_row), "set_row: row out of bounds");
    
    new_row = Matrix<double>({
        {1, 2, 3, 4, 5, 6},
        {1, 2, 3, 4, 5, 6}
    });
    XCTAssertThrows(result.set_row(1, new_row), "set_row: new row has more than one row");
    
    new_row = Matrix<double>({
        {1, 2, 3, 4, 5, 6}
    });
    XCTAssertThrows(result.set_row(1, new_row), "set_row: new row has more elements than old one");
    
    result.ones();
    new_row = Matrix<double>({{2, -5, 7, 9, -7}});
    result.set_row(1, new_row.get_row_as_vector(0));
    XCTAssertEqual(result, expected, "set_row (vector) complete OK");
    
    XCTAssertThrows(result.set_row(4, new_row.get_row_as_vector(0)), "set_row (vector): row out of bounds");
    
    new_row = Matrix<double>({{1, 2, 3, 4, 5, 6}});
    XCTAssertThrows(result.set_row(1, new_row.get_row_as_vector(0)), "set_row (vector): new row has more elements than old one");
    
    new_row = Matrix<double>({{1, 2, 3}});
    XCTAssertThrows(result.set_row(1, new_row.get_row_as_vector(0)), "set_row (vector): new row has fewer elements than old one");
}

- (void)testSetColumn {
    auto result = Matrix<double>::ones(4, 5);
    
    const Matrix<double> expected({
        { 1,  2, 1, 1, 1},
        { 1, -3, 1, 1, 1},
        { 1,  6, 1, 1, 1},
        { 1, -8, 1, 1, 1}
    });
    
    Matrix<double> new_column(4, 1, {2, -3, 6, -8});
    result.set_column(1, new_column);
    XCTAssertEqual(result, expected, "set_column complete OK");
    
    XCTAssertThrows(result.set_column(5, new_column), "set_column: column out of bounds");
    
    new_column = Matrix<double>(4, 2, {1, 2, 3, 4, 5, 6, 7, 8});
    XCTAssertThrows(result.set_column(1, new_column), "set_column: new column has more than one column");
    
    new_column = Matrix<double>(5, 1, {1, 2, 3, 4, 5});
    XCTAssertThrows(result.set_column(1, new_column), "set_column: new column has more elements than old one");
    
    result.ones();
    new_column = Matrix<double>(4, 1, {2, -3, 6, -8});
    result.set_column(1, new_column.get_column_as_vector(0));
    XCTAssertEqual(result, expected, "set_column (vector) complete OK");
    
    XCTAssertThrows(result.set_column(5, new_column.get_column_as_vector(0)), "set_column (vector): column out of bounds");
    
    new_column = Matrix<double>(5, 1, {1, 2, 3, 4, 5});
    XCTAssertThrows(result.set_column(1, new_column.get_column_as_vector(0)), "set_column (vector): new column has more elements than old one");
    
    new_column = Matrix<double>(3, 1, {1, 2, 3});
    XCTAssertThrows(result.set_column(1, new_column.get_column_as_vector(0)), "set_column (vector): new column has fewer elements than old one");
}

- (void)testSetDiagonal {
    auto result1 = Matrix<double>::ones(4, 4);
    result1.SetDiagonal(2);
    
    const Matrix<double> expected1({
        {2, 1, 1, 1},
        {1, 2, 1, 1},
        {1, 1, 2, 1},
        {1, 1, 1, 2}
    });
    
    XCTAssertEqual(result1, expected1, "SetDiagonal scalar value OK");
    
    auto result2 = Matrix<double>::zero(5, 5);
    result2.SetDiagonal(Vector<double>({1, 2, 3, 4, 5}));
    
    const Matrix<double> expected2({
        {1, 0, 0, 0, 0},
        {0, 2, 0, 0, 0},
        {0, 0, 3, 0, 0},
        {0, 0, 0, 4, 0},
        {0, 0, 0, 0, 5}
    });
    
    XCTAssertEqual(result2, expected2, "SetDiagonal vector values OK");
    
    Matrix<double> result3(4, 5);
    XCTAssertThrows(result3.SetDiagonal(4), "Matrix must be square");
    XCTAssertThrows(result3.SetDiagonal(Vector<double>(4)), "Matrix must be square");
    XCTAssertThrows(result1.SetDiagonal(Vector<double>(5)), "The number of rows of the matrix and the size of the vector does not match");
}

- (void)testSumRows {
    Matrix<double> matrix_test({
        { 3,  2,  1,  2},
        { 7,  6,  5,  1},
        {12, 10,  9,  8},
        {15, 14, 13, 12}
    });
    
    Matrix<double> expected1(4, 1, {8, 19, 39, 54});
    XCTAssertEqual(matrix_test.SumRows(), expected1, "SumRows OK");
    XCTAssertEqual(matrix_test.SumRowsAsVector(), expected1.get_column_as_vector(0), "SumRowsAsVector OK");
    
    
    XCTAssertEqual(matrix_test.SumRow(1), expected1[0][1], "SumRow OK");
    XCTAssertThrows(matrix_test.SumRow(4), "SumRow out of range");
}

- (void)testSumColumns {
    const Matrix<double> matrix_test({
        { 3,  2,  1,  2},
        { 7,  6,  5,  1},
        {12, 10,  9,  8},
        {15, 14, 13, 12}
    });
    
    const Matrix<double> expected1(1, 4, {37, 32, 28, 23});
    XCTAssertEqual(matrix_test.SumColumns(), expected1, "CumColumns OK");
    XCTAssertEqual(matrix_test.SumColumnsAsVector(), expected1.get_row_as_vector(0), "SumColumnsAsVector OK");
    
    XCTAssertEqual(matrix_test.SumColumn(1), expected1[0][1], "SumColumn OK");
    XCTAssertThrows(matrix_test.SumColumn(4), "SumColumn out of range");
}

- (void)testIfstreamOperator {
    
    std::ifstream file("data/math/containers/BigMatrix.csv", std::ios::in);
    
    // Big Matrix
    Matrix<double> matrix;
    file >> matrix;
    file.close();
    
    std::pair<size_t, size_t> expected_dimensions {449, 106};
    
    XCTAssertEqual(matrix.Dimensions(), expected_dimensions, "Big Matrix dimensions are the right ones.");
    
    // Small Matrix
    file.open("data/math/containers/SmallMatrix.csv", std::ios::in);
    file >> matrix;
    file.close();
    
    Matrix<double> expected_matrix({
        {0.1680625202081561,  0.1722621842917498,  0.7412169766918424,  0.6987185197938875,  0.3302779764663697},
        {0.10215466928767196, 0.3990300364707511,  0.7262335341926227,  0.08155865621143804, 0.3684479022981741},
        {0.4550947968717964,  0.33873967383847237, 0.4988407455385848,  0.8256508889575926,  0.4998906510004011},
        {0.6474657582972254,  0.5223187955808917,  0.548139118672313,   0.7215750817149178,  0.14924930831433234},
        {0.796918578539592,   0.4904564565638245,  0.11938391918190965, 0.9765232400263497,  0.6245631592365628}
    });
    
    XCTAssertEqual(matrix, expected_matrix, "Small Matrix has been loaded properly.");
    
    std::ifstream unavailable_file("data/unavailable_file", std::ios::in);
    Matrix<double> matrix1;
    
    XCTAssertThrows(unavailable_file >> matrix1, "File does not exist OK");
    unavailable_file.close();
}

- (void)testProductBetweenMatrixAndVector {
    Matrix<double> matrix({
        { 3},
        { 7},
        {12},
        {15}
    });
    
    Vector<double> vector({ 7,  6,  5,  1});
    
    const Matrix<double> expected({
        { 21, 18, 15,  3},
        { 49, 42, 35,  7},
        { 84, 72, 60, 12},
        {105, 90, 75, 15}
    });
    
    XCTAssertEqual(matrix * vector, expected, "Product between matrix and vector OK");
    
    Matrix<double> matrix_2columns({
        { 3, 4},
        { 7, 1},
        {12, 5},
        {15, 6}
    });
    
    XCTAssertThrows(matrix_2columns * vector, "Dimensions are not compatible");
}

- (void)testProductBetweenVectorAndMatrix {
    const Vector<double> vector({ 7,  6,  5,  1});
    const Matrix<double> matrix({
        { 3, 4},
        { 7, 8},
        {12, 20},
        {15, 3}
    });
    
    XCTAssertEqual(vector * matrix, Vector<double>({138, 179}), "Product between vector and matrix OK");
    
    XCTAssertThrows(Vector<double>({1, 35, 362, 4362, 34}) * matrix, "Dimensions are not compatible");
}

- (void)testTransposeVectorToMatrix {
    const Vector<double> vector({0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
    const Matrix<double> expected({
        {0},
        {1},
        {2},
        {3},
        {4},
        {5},
        {6},
        {7},
        {8},
        {9}
    });
    
    XCTAssertEqual(Transpose(vector), expected, "Transpose vector is the expected matrix");
}

- (void)testMaximumAndMinimumElements {
    const Matrix<double> matrix({
        {  1.134,   0.001,   2.523,  -0.231,     0.321, -312353.123},
        {  5.213,   6.312,  -7.142,   8.243,     9.234,      21.426},
        { 10.123, -11.321,  12.213, -13.213,    14.231,       1.213},
        { -0.024,   0.314,  17.143,  18.143, -1913.136,  523251.316},
        {  5.432,  -6.236,   7.342,   8.324,    -9.341,      21.341}
    });
    
    XCTAssertEqual(matrix.max_element(), 523251.316, "Maximum element OK");
    XCTAssertEqual(matrix.abs_max_element(), 523251.316, "Absolute maximum element OK");
    XCTAssertEqual(matrix.abs_max_element_with_sign(), 523251.316, "Absolute maximum element with sign OK");
    
    XCTAssertEqual(matrix.min_element(), -312353.123, "Minimum element OK");
    XCTAssertEqual(matrix.abs_min_element(), 0.001, "Absolute minimum element OK");
    XCTAssertEqual(matrix.abs_min_element_with_sign(), 0.001, "Absolute minimum element with sign OK");
}
    
- (void)testFindMethod {
    const Matrix<double> matrix({{1, 2, 3}, {4, 5, 6}, {7, 4, 9}});
    
    auto it_value_4 = matrix.find(4);
    XCTAssertEqual(*it_value_4, 4, "Value 4 has been found");
    
    auto it_first_value_4 = matrix.begin() + 3;
    XCTAssertEqual(it_value_4, it_first_value_4, "Found value 4 is the firstone");
    
    auto it_value_12 = matrix.find(12);
    XCTAssertEqual(it_value_12, matrix.end(), "Value 12 is not present");
}

- (void)testAccessors {
    const Matrix<double> matrix1({
        { 21, 18, 15,  3},
        { 49, 42, 35,  7},
        { 84, 72, 60, 12},
        {105, 90, 75, 15}
    });
    
    XCTAssertEqual(matrix1.at(0, 0), 21, "Matrix at for first element OK");
    XCTAssertEqual(matrix1.at(3, 3), 15, "Matrix at for last element OK");
    XCTAssertEqual(matrix1.at(0, 3), 3, "Matrix at for last element of first row OK");
    XCTAssertEqual(matrix1.at(3, 0), 105, "Matrix at for first element of last row OK");
    
    XCTAssertThrows(matrix1.at(-1, -1), "Matrix at out of range for negative indexes");
    XCTAssertThrows(matrix1.at(4, 4), "Matrix at out of range for indexes greater than matrix dimensions");
    XCTAssertThrows(matrix1.at(0, 4), "Matrix at out of range for index that exceeds the limits of the row");
    XCTAssertThrows(matrix1.at(4, 0), "Matrix at out of range for index that exceeds the limits of the column");
    
    Matrix<double> matrix2({
        { 21, 18, 15,  3},
        { 49, 42, 35,  7},
        { 84, 72, 60, 12},
        {105, 90, 75, 15}
    });
    
    matrix2.at(0, 0) = 12;
    XCTAssertEqual(matrix2.at(0, 0), 12, "Matrix at for setting the first element OK");
    
    matrix2.at(3, 3) = 51;
    XCTAssertEqual(matrix2.at(3, 3), 51, "Matrix at for setting the last element OK");
    
    matrix2.at(0, 3) = -3;
    XCTAssertEqual(matrix2.at(0, 3), -3, "Matrix at for setting the last element of first row OK");
    
    matrix2.at(3, 0) = 501;
    XCTAssertEqual(matrix2.at(3, 0), 501, "Matrix at for setting the first element of last row OK");
    
    XCTAssertThrows(matrix2.at(-1, -1) = 21, "Matrix at out of range for setting value at negative indexes");
    XCTAssertThrows(matrix2.at(4, 4) = 15, "Matrix at out of range for setting value at indexes greater than matrix dimensions");
    XCTAssertThrows(matrix2.at(0, 4) = 3, "Matrix at out of range for setting value at index that exceeds the limits of the row");
    XCTAssertThrows(matrix2.at(4, 0) = 105, "Matrix at out of range for setting value at index that exceeds the limits of the column");
}

- (void)testIsNullAndIsEmptyMethods {
    Matrix<double> matrix;
    XCTAssert(matrix.is_null(), "Matrix is null");
    XCTAssert(matrix.is_empty(), "Matrix is empty");
    
    matrix = Matrix<double>::zero(4, 4);
    XCTAssert(matrix.is_null(), "Matrix is null");
    XCTAssert(!matrix.is_empty(), "Matrix is not empty");
}

- (void)testClear {
    Matrix<double> matrix({
        { 21, 18, 15,  3},
        { 49, 42, 35,  7},
        { 84, 72, 60, 12},
        {105, 90, 75, 15}
    });
    
    matrix.clear();
    
    XCTAssert(matrix.is_null(), "Matrix is null");
    XCTAssert(matrix.is_empty(), "Matrix is empty");
}

- (void)testHasDuplicate {
    Matrix<double> matrix({
        { 21, 18, 15,  3},
        { 49, 42, 35,  7},
        { 84, 72, 60, 12},
        {105, 90, 75, 15}
    });
    
    XCTAssert(matrix.has_duplicate(), "Matrix has duplicate elements");
    
    matrix = Matrix<double>({
        { 21, 18, 15,  3},
        { 49, 42, 35,  7},
        { 84, 72, 60, 12},
        {105, 90, 75, 16}
    });
    
    XCTAssert(!matrix.has_duplicate(), "Matrix does not have duplicate elements");
    
    Matrix<double> matrix2({
        { 21.1727, 18.1355, 15.5478,  3.5235},
        { 49.1253, 42.2135, 35.2135,  7.2135},
        { 84.3255, 72.3164, 60.2544, 12.7543},
        {105.2757, 90.3164, 75.7235, 15.5435}
    });
    
    XCTAssert(!matrix2.has_duplicate(), "Matrix double does not have duplicate elements with accuracy inf");
    XCTAssert(!matrix2.has_duplicate(1E-04), "Matrix double does not have duplicate elements with accuracy 1E-04");
    XCTAssert(matrix2.has_duplicate(1E-02), "Matrix double does not have duplicate elements with accuracy 1E-02");
}

- (void)testDeterminant {
    const Matrix<double> matrix1({
        { 21, 18, 15,  4},
        { 49, 41, 35,  7},
        { 84, 72, 63, 12},
        {105, 90, 75, 15}
    });
    
    XCTAssertEqual(matrix1.Determinant(), 315, "Determinant OK");
    
    const Matrix<double> matrix2({
        { 21, 18, 15,  4},
        { 42, 36, 30,  8},
        { 84, 72, 63, 12},
        {105, 90, 75, 15}
    });
    
    XCTAssertEqual(matrix2.Determinant(), 0, "Determinant 0 OK");
    
    const Matrix<double> matrix3({
        { 21, 18, 15,  4},
        { 42, 36, 30,  8},
        { 84, 72, 63, 12}
    });
    
    XCTAssertThrows(matrix3.Determinant(), "Unable to compute determinant of a non-square matrix");
}

- (void)testOfstreamOperator {
    
    std::string expected_content;
    {
        std::ifstream expected_file;
        expected_file.open("data/math/containers/MatrixTests-testOfstreamOperator.txt", std::ios::in);
        expected_content = std::string(std::istreambuf_iterator<char>(expected_file),
                                       std::istreambuf_iterator<char>());
        expected_file.close();
    }
    
    // --- Test
    std::stringstream test_output;
    std::streambuf *coutbuf = std::cout.rdbuf();
    
    // Redirect std::cout buffer
    std::cout.rdbuf(test_output.rdbuf());
    
    const Matrix<double> matrix1({
        {21, 18, 15,  4},
        {42, 36, 30,  8},
        {84, 72, 63, 12}
    });
    std::cout << matrix1;
    
    const Matrix<double> matrix2({
        {21, 18, 15,  4}
    });
    std::cout << matrix2;
    
    const Matrix<double> matrix3({
        {1}
    });
    std::cout << matrix3;
    
    // Restore std::cout buffer
    std::cout.rdbuf(coutbuf);
    test_output << matrix1;
    
    XCTAssertEqual(test_output.str(), expected_content, "Ofstream operator OK");
}

@end
