//
//  MatrixTests.m
//  Tests
//
//  Created by Carlos David on 11/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#import <XCTest/XCTest.h>

#import "../../TestsTools.h"
#import "../../../NumericalPDEs/math/containers/matrix.hpp"

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
    
    const Matrix<int> matrix;
    XCTAssert(matrix.IsEmpty(), "matrix is empty");
    XCTAssert(matrix.IsNull(), "matrix is null");
    
    const Matrix<int> matrix1(4, 4, 1);
    XCTAssert(!matrix1.IsEmpty() && !matrix1.IsNull(), "matrix1 is not empty and is not null");
    XCTAssertEqual(matrix1.Rows(), 4, "The number of rows of matrix1 is OK");
    XCTAssertEqual(matrix1.Columns(), 4, "The number of columns of matrix1 is OK");
    XCTAssertEqual(matrix1.Size(), 16, "The size of matrix1 is OK");
    
    bool all_elements_are_one = true;
    for (auto it = matrix1.Begin(); it != matrix1.End(); ++it) {
        if (*it != 1) {
            all_elements_are_one = false;
            break;
        }
    }
    XCTAssert(all_elements_are_one, "All elements of matrix1 are OK");
    
    const Matrix<int> matrix2(4, 4, {
        0,  1,  2,  3,
        4,  5,  6,  7,
        8,  9, 10, 11,
        12, 13, 14, 15
    });
    XCTAssert(!matrix2.IsEmpty() && !matrix2.IsNull(), "matrix2 is not empty and is not null");
    XCTAssertEqual(matrix2.Rows(), 4, "The number of rows of matrix2 is OK");
    XCTAssertEqual(matrix2.Columns(), 4, "The number of columns of matrix2 is OK");
    XCTAssertEqual(matrix2.Size(), 16, "The size of matrix2 is OK");
    
    XCTAssertThrows(Matrix<int>(2, 2, { 0,  1,  2 }), "The array of values has different size than the matrix");
    
    const Matrix<int> matrix3({
        { 0,  1,  2,  3},
        { 4,  5,  6,  7},
        { 8,  9, 10, 11},
        {12, 13, 14, 15}
    });
    XCTAssert(!matrix3.IsEmpty() && !matrix3.IsNull(), "matrix3 is not empty and is not null");
    XCTAssertEqual(matrix3.Rows(), 4, "The number of rows of matrix3 is OK");
    XCTAssertEqual(matrix3.Columns(), 4, "The number of columns of matrix3 is OK");
    XCTAssertEqual(matrix3.Size(), 16, "The size of matrix3 is OK");
    
    Matrix<int> matrix4(matrix3);
    XCTAssertEqual(matrix4, matrix3, "matrix3 and matrix4 are equal");
    
    Matrix<int> matrix5(std::move(matrix4));
    XCTAssertEqual(matrix5, matrix3, "matrix3 and matrix4 are equal");
    XCTAssert(matrix4.IsEmpty(), "matrix4 is empty after had been moved");
    XCTAssert(matrix4.IsNull(), "matrix4 is null after had been moved");
    
    matrix4 = matrix2;
    XCTAssertEqual(matrix4, matrix2, "matrix4 and matrix2 are equal");
    
    matrix5 = std::move(matrix4);
    XCTAssertEqual(matrix5, matrix2, "matrix5 and matrix2 are equal");
    XCTAssert(matrix4.IsEmpty(), "matrix4 is empty after had been moved in asignation");
    XCTAssert(matrix4.IsNull(), "matrix4 is null after had been moved in asignation");
}

- (void)testOnes {
    const Matrix<uint> expected_ones(5, 5, 1);
    XCTAssertEqual(Matrix<uint>::Ones(5, 5), expected_ones, "Ones static method OK");
    
    Matrix<uint> result(5, 5, 2);
    result.Ones();
    XCTAssertEqual(result, expected_ones, "Ones method OK");
}

- (void)testZero {
    const Matrix<uint> expected_zero(5, 5, 0);
    XCTAssertEqual(Matrix<uint>::Zero(5, 5), expected_zero, "Zero static method OK");
    
    Matrix<uint> result(5, 5, 2);
    result.Zero();
    XCTAssertEqual(result, expected_zero, "Zero method OK");
}

- (void)testIdentity {
    const Matrix<uint> expected_identity({
        {1, 0, 0, 0, 0},
        {0, 1, 0, 0, 0},
        {0, 0, 1, 0, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 0, 0, 1}
    });
    XCTAssertEqual(Matrix<uint>::Identity(5, 5), expected_identity, "Identity static method OK");
    
    Matrix<uint> result(5, 5, 2);
    result.Identity();
    XCTAssertEqual(result, expected_identity, "Identity method OK");
    
    Matrix<uint> result2(5, 4, 2);
    XCTAssertThrows(result2.Identity(), "Identity method is only valid for squared matrices");
}

- (void)testComparison {
    const Matrix<int> matrix1({
        { 0,  1,  2,  3},
        { 4,  5,  6,  7},
        { 8,  9, 10, 11},
        {12, 13, 14, 15}
    });
    
    XCTAssert(matrix1 == matrix1, "Matrix is equal to itself");
    XCTAssert(!(matrix1 != matrix1), "Matrix is not different to itself");
    
    const Matrix<int> matrix2({
        { 0,  1,  2,  3},
        { 8,  9, 10, 11},
        { 4,  5,  6,  7},
        {12, 13, 14, 15}
    });
    
    XCTAssert(matrix1 != matrix2, "matrix1 is different to matrix2");
    XCTAssert(!(matrix1 == matrix2), "matrix1 is not equal to matrix2");
    
    // Same as matrix1 whit less columns
    const Matrix<int> matrix3({
        { 0,  1,  2},
        { 4,  5,  6},
        { 8,  9, 10},
        {12, 13, 14}
    });
    
    XCTAssert(matrix3 != matrix1, "matrix3 is different to matrix1");
    XCTAssert(!(matrix3 == matrix1), "matrix3 is not equal to matrix1");
    
    // Same as matrix2 whit less rows
    const Matrix<int> matrix4({
        { 0,  1,  2,  3},
        { 4,  5,  6,  7},
        {12, 13, 14, 15}
    });
    
    XCTAssert(matrix4 != matrix2, "matrix4 is different to matrix2");
    XCTAssert(!(matrix4 == matrix2), "matrix4 is not equal to matrix2");
}

- (void)testResize {
    const Matrix<int> matrix({
        {0, 1,  2},
        {4, 5,  6},
        {8, 9, 10}
    });
    
    auto result(matrix);
    result.Resize(matrix.Rows(), matrix.Columns());
    XCTAssertEqual(result, matrix, "matrix has not been changed");
    
    const Matrix<int> expected1({
        {0, 1},
        {4, 5}
    });
    
    result = Matrix<int>(matrix);
    result.Resize(2, 2);
    XCTAssertEqual(result, expected1, "matrix has been resized OK to a 2x2 matrix");
    
    const Matrix<int> expected2({
        {0, 1},
        {4, 5},
        {8, 9}
    });
    
    result = Matrix<int>(matrix);
    result.Resize(3, 2);
    XCTAssertEqual(result, expected2, "matrix has been resized OK to a 3x2 matrix");
    
    const Matrix<int> expected3({
        {0, 1, 2},
        {4, 5, 6}
    });
    
    result = Matrix<int>(matrix);
    result.Resize(2, 3);
    XCTAssertEqual(result, expected3, "matrix has been resized OK to a 2x3 matrix");
    
    const Matrix<int> expected4({
        {0, 1,  2, 0},
        {4, 5,  6, 0},
        {8, 9, 10, 0},
        {0, 0,  0, 0}
    });
    
    result = Matrix<int>(matrix);
    result.Resize(4, 4, true);
    XCTAssertEqual(result, expected4, "matrix has been resized OK to a 4x4 matrix adding zeros to the new elements");
    
    const Matrix<int> expected5({
        {0, 1,  2},
        {4, 5,  6},
        {8, 9, 10},
        {0, 0,  0}
    });
    
    result = Matrix<int>(matrix);
    result.Resize(4, 3, true);
    XCTAssertEqual(result, expected5, "matrix has been resized OK to a 4x3 matrix adding zeros to the new elements");
}

- (void)testChangeDimensions {
    const Matrix<int> matrix({
        {0, 1,  2, 4},
        {4, 5,  6, 7},
        {8, 9, 10, 11}
    });
    
    const Matrix<int> expected({
        {0,  1,  2},
        {4,  4,  5},
        {6,  7,  8},
        {9, 10, 11}
    });
    
    auto result(matrix);
    result.ChangeDimensions(4, 3);
    
    XCTAssertEqual(result, expected, "Dimensions have been changed OK");
    XCTAssertThrows(result.ChangeDimensions(3, 3), "The total size of the matrix cannot be changed");
}

- (void)testProducts {
    const Matrix<int> matrix1({
        { 0,  1,  2,  3},
        { 4,  5,  6,  7},
        { 8,  9, 10, 11},
        {12, 13, 14, 15}
    });
    
    const Matrix<int> matrix2({
        { 3,  2,  1,  0},
        { 7,  6,  5,  4},
        {11, 10,  9,  8},
        {15, 14, 13, 12}
    });
    
    const Matrix<int> expected1({
        { 74,  68,  62,  56},
        {218, 196, 174, 152},
        {362, 324, 286, 248},
        {506, 452, 398, 344}
    });
    
    XCTAssertEqual(matrix1 * matrix2, expected1, "Square matrices product OK");
    
    const Matrix<int> matrix3({
        { 0,  1,  2,  3, -1},
        { 4,  5,  6,  7, -2},
        { 8,  9, 10, 11, -3},
        {12, 13, 14, 15, -4}
    });
    
    const Matrix<int> matrix4({
        { 3,   2,   1,   0, -11, 235},
        { 7,   6,   5,   4, -22, 264},
        {11,  10,   9,   8, -33, 436},
        {15,  14,  13,  12, -44, 643},
        {54, 235,  21,  21, 235, 426}
    });
    
    const Matrix<int> expected2({
        { 20, -167,  41,  35,  -455,  2639},
        {110, -274, 132, 110, -1130,  8525},
        {200, -381, 223, 185, -1805, 14411},
        {290, -488, 314, 260, -2480, 20297}
    });
    
    XCTAssertEqual(matrix3 * matrix4, expected2, "Rectangular matrices product OK");
    
    XCTAssertThrows(matrix4 * matrix3, "Matrices dimensions are incompatible");
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
    
    XCTAssertEqual(matrix.Transpose(), expected, "Transpose operation OK");
}

- (void)testSquaredMatrix {
    const auto squared_matrix = Matrix<int>::Ones(4, 4);
    XCTAssert(squared_matrix.IsSquared(), "squared_matrix is squared");
    
    const auto rectangular_matrix = Matrix<int>::Ones(4, 3);
    XCTAssert(!rectangular_matrix.IsSquared(), "rectangular_matrix is not squared");
}

- (void)testGettersMethods {
    const Matrix<int> matrix({
        {0,  1,  2,  3,  4},
        {5,  6,  7,  8,  9},
        {10, 11, 12, 13, 14},
        {15, 16, 17, 18, 19}
    });
    
    // --- GetMatrix
    const Matrix<int> expected_matrix1({
        {7,  8,  9},
        {12, 13, 14},
        {17, 18, 19}
    });
    
    XCTAssertEqual(matrix.GetMatrix(1, 2), expected_matrix1, "GetMatrix without lengths OK");
    
    const Matrix<int> expected_matrix2({
        {0,  1,  2},
        {5,  6,  7},
        {10, 11, 12}
    });
    
    XCTAssertEqual(matrix.GetMatrix(0, 0, 3, 3), expected_matrix2, "GetMatrix with lengths OK");
    
    XCTAssertThrows(matrix.GetMatrix(1, 1, 6, 6), "Unable to get submatrix. Elements out of bounds");
    
    // --- GetRow
    const Matrix<int> expected_first_row({{0,  1,  2,  3,  4}});
    XCTAssertEqual(matrix.GetRow(0), expected_first_row, "GetRow first OK");
    
    const Matrix<int> expected_last_row({{15, 16, 17, 18, 19}});
    XCTAssertEqual(matrix.GetRow(3), expected_last_row, "GetRow last OK");
    
    const Matrix<int> expected_in_between_row({{10, 11, 12, 13, 14}});
    XCTAssertEqual(matrix.GetRow(2), expected_in_between_row, "GetRow last OK");
    
    XCTAssertThrows(matrix.GetRow(4), "GetRow out of bounds");
    
    // --- GetColumn
    const Matrix<int> expected_first_column(4, 1, {0, 5, 10, 15});
    XCTAssertEqual(matrix.GetColumn(0), expected_first_column, "GetRow first OK");
    
    const Matrix<int> expected_last_column(4, 1, {4, 9, 14, 19});
    XCTAssertEqual(matrix.GetColumn(4), expected_last_column, "GetRow last OK");
    
    const Matrix<int> expected_in_between_column(4, 1, {2, 7, 12, 17});
    XCTAssertEqual(matrix.GetColumn(2), expected_in_between_column, "GetRow last OK");
    
    XCTAssertThrows(matrix.GetColumn(5), "GetColumn out of bounds");
    
    // --- GetDiagonal
    XCTAssertThrows(matrix.GetDiagonal(), "This method is only available for squared matrices");
    
    const Matrix<int> squared_matrix({
        { 0,  1,  2,  3},
        { 5,  6,  7,  8},
        {10, 11, 12, 13},
        {15, 16, 17, 18}
    });
    
    const Vector<int> expected_diagonal({0, 6, 12, 18});
    XCTAssertEqual(squared_matrix.GetDiagonal(), expected_diagonal, "GetDiagonal method OK");
}

- (void)testGettersAsVectorMethods {
    const Matrix<int> matrix({
        {0,  1,  2,  3,  4},
        {5,  6,  7,  8,  9},
        {10, 11, 12, 13, 14},
        {15, 16, 17, 18, 19}
    });
    
    // --- GetColumnAsVector ---
    // Get the whole first column
    Vector<int> expected_vector({0, 5, 10, 15});
    auto result_vector = matrix.GetColumnAsVector(0);
    XCTAssertEqual(result_vector, expected_vector, "GetColumnAsVector for the whole first column OK!");
    
    // Get the last two elements of the first column
    expected_vector = Vector<int>({10, 15});
    result_vector = matrix.GetColumnAsVector(0, 2);
    XCTAssertEqual(result_vector, expected_vector, "GetColumnAsVector for the last two elements of the first column OK!");
    
    // Get the whole third column
    expected_vector = Vector<int>({2, 7, 12, 17});
    result_vector = matrix.GetColumnAsVector(2);
    XCTAssertEqual(result_vector, expected_vector, "GetColumnAsVector for the whole third column OK!");
    
    // Get the last three elements of the fourth column
    expected_vector = Vector<int>({8, 13, 18});
    result_vector = matrix.GetColumnAsVector(3, 1);
    XCTAssertEqual(result_vector, expected_vector, "GetColumnAsVector for the last three elements of the fourth column OK!");
    
    // Get the last column
    expected_vector = Vector<int>({4, 9, 14, 19});
    result_vector = matrix.GetColumnAsVector(4);
    XCTAssertEqual(result_vector, expected_vector, "GetColumnAsVector for the last column OK!");
    
    // Get column out of bounds
    XCTAssertThrows(matrix.GetColumnAsVector(3, 6), "GetColumnAsVector out of bounds by number of elements");
    XCTAssertThrows(matrix.GetColumnAsVector(6), "GetColumnAsVector out of bounds by number of column");
    
    // --- GetRowAsVector ---
    // Get the whole first row
    expected_vector = Vector<int>({0,  1,  2,  3,  4});
    result_vector = matrix.GetRowAsVector(0);
    XCTAssertEqual(result_vector, expected_vector, "GetRowAsVector for the whole first row OK!");
    
    // Get the last two elements of the first column
    expected_vector = Vector<int>({7,  8,  9});
    result_vector = matrix.GetRowAsVector(1, 2);
    XCTAssertEqual(result_vector, expected_vector, "GetRowAsVector for the last three elements of the second row OK!");
    
    // Get the whole third row
    expected_vector = Vector<int>({10, 11, 12, 13, 14});
    result_vector = matrix.GetRowAsVector(2);
    XCTAssertEqual(result_vector, expected_vector, "GetRowAsVector for the whole third row OK!");
    
    // Get the last element of the last row
    expected_vector = Vector<int>({19});
    result_vector = matrix.GetRowAsVector(3, 4);
    XCTAssertEqual(result_vector, expected_vector, "GetRowAsVector for the last element of the last row OK!");
    
    // Get row out of bounds
    XCTAssertThrows(matrix.GetRowAsVector(3, 6), "GetColumnAsVector out of bounds by number of elements");
    XCTAssertThrows(matrix.GetRowAsVector(6), "GetColumnAsVector out of bounds by number of column");
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
    
    XCTAssertEqual(result, expected1, "SetMatrix without lengths OK");
    
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
    
    XCTAssertEqual(result, expected2, "SetMatrix with lengths OK");
    
    result.SetMatrix(0, 0,
                     Matrix<int>({
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
    auto result = Matrix<int>::Ones(4, 5);
    
    const Matrix<int> expected({
        { 1,  1, 1, 1,  1},
        { 2, -5, 7, 9, -7},
        { 1,  1, 1, 1,  1},
        { 1,  1, 1, 1,  1}
    });
    
    Matrix<int> new_row({{2, -5, 7, 9, -7}});
    result.SetRow(1, new_row);
    XCTAssertEqual(result, expected, "SetRow complete OK");
    
    XCTAssertThrows(result.SetRow(4, new_row), "SetRow: row out of bounds");
    
    new_row = Matrix<int>({
        {1, 2, 3, 4, 5, 6},
        {1, 2, 3, 4, 5, 6}
    });
    XCTAssertThrows(result.SetRow(1, new_row), "SetRow: new row has more than one row");
    
    new_row = Matrix<int>({
        {1, 2, 3, 4, 5, 6}
    });
    XCTAssertThrows(result.SetRow(1, new_row), "SetRow: new row has more elements than old one");
    
    result.Ones();
    new_row = Matrix<int>({{2, -5, 7, 9, -7}});
    result.SetRow(1, new_row.GetRowAsVector(0));
    XCTAssertEqual(result, expected, "SetRow (vector) complete OK");
    
    XCTAssertThrows(result.SetRow(4, new_row.GetRowAsVector(0)), "SetRow (vector): row out of bounds");
    
    new_row = Matrix<int>({{1, 2, 3, 4, 5, 6}});
    XCTAssertThrows(result.SetRow(1, new_row.GetRowAsVector(0)), "SetRow (vector): new row has more elements than old one");
    
    new_row = Matrix<int>({{1, 2, 3}});
    XCTAssertThrows(result.SetRow(1, new_row.GetRowAsVector(0)), "SetRow (vector): new row has fewer elements than old one");
}

- (void)testSetColumn {
    auto result = Matrix<int>::Ones(4, 5);
    
    const Matrix<int> expected({
        { 1,  2, 1, 1, 1},
        { 1, -3, 1, 1, 1},
        { 1,  6, 1, 1, 1},
        { 1, -8, 1, 1, 1}
    });
    
    Matrix<int> new_column(4, 1, {2, -3, 6, -8});
    result.SetColumn(1, new_column);
    XCTAssertEqual(result, expected, "SetColumn complete OK");
    
    XCTAssertThrows(result.SetColumn(5, new_column), "SetColumn: column out of bounds");
    
    new_column = Matrix<int>(4, 2, {1, 2, 3, 4, 5, 6, 7, 8});
    XCTAssertThrows(result.SetColumn(1, new_column), "SetColumn: new column has more than one column");
    
    new_column = Matrix<int>(5, 1, {1, 2, 3, 4, 5});
    XCTAssertThrows(result.SetColumn(1, new_column), "SetColumn: new column has more elements than old one");
    
    result.Ones();
    new_column = Matrix<int>(4, 1, {2, -3, 6, -8});
    result.SetColumn(1, new_column.GetColumnAsVector(0));
    XCTAssertEqual(result, expected, "SetColumn (vector) complete OK");
    
    XCTAssertThrows(result.SetColumn(5, new_column.GetColumnAsVector(0)), "SetColumn (vector): column out of bounds");
    
    new_column = Matrix<int>(5, 1, {1, 2, 3, 4, 5});
    XCTAssertThrows(result.SetColumn(1, new_column.GetColumnAsVector(0)), "SetColumn (vector): new column has more elements than old one");
    
    new_column = Matrix<int>(3, 1, {1, 2, 3});
    XCTAssertThrows(result.SetColumn(1, new_column.GetColumnAsVector(0)), "SetColumn (vector): new column has fewer elements than old one");
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
    const Matrix<int> matrix_test({
        { 3,  2,  1,  2},
        { 7,  6,  5,  1},
        {12, 10,  9,  8},
        {15, 14, 13, 12}
    });
    
    const Matrix<int> expected1(1, 4, {37, 32, 28, 23});
    auto result1 = matrix_test.SumColumns();
    
    XCTAssert(result1 == expected1, "CumColumns OK");
    
    auto expected2 = expected1.GetRowAsVector(0);
    auto result2 = matrix_test.SumColumnsAsVector();
    
    XCTAssert(result2 == expected2, "SumColumnsAsVector OK");
}

- (void)testLoadMatrixFromFile {
    
    std::ifstream file;
    Matrix<double> matrix;
    
    // Big Matrix
    file.open("data/math/containers/BigMatrix.csv", std::ios::in);
    file >> matrix;
    file.close();
    
    std::pair<size_t, size_t> expected_dimensions {449, 106};
    
    XCTAssert(matrix.Dimensions() == expected_dimensions, "Big Matrix dimensions are the right ones.");
    
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
    
    XCTAssert(matrix == expected_matrix, "Small Matrix has been loaded properly.");
}

- (void)testProduct {
    const Matrix<int> matrix1({
        { 3,  2,  1,  2},
        { 7,  6,  5,  1},
        {12, 10,  9,  8},
        {15, 14, 13, 12}
    });
    
    const Matrix<int> matrix2({
        { 2,  1,  2},
        { 6,  5,  1},
        {10,  9,  8},
        {14, 13, 12}
    });
    
    const Matrix<int> expected({
        { 56,  48,  40},
        {114,  95,  72},
        {286, 247, 202},
        {412, 358, 292}
    });
    
    auto result = matrix1 * matrix2;
    
    XCTAssertEqual(result, expected, "Product of two patrices OK");
    XCTAssertThrows(matrix2 * matrix1, "Incompatible dimensions to compute the product");
}

- (void)testProductBetweenMatrixAndVector {
    Matrix<int> matrix({
        { 3},
        { 7},
        {12},
        {15}
    });
    
    Vector<int> vector({ 7,  6,  5,  1});
    
    const Matrix<int> expected({
        { 21, 18, 15,  3},
        { 49, 42, 35,  7},
        { 84, 72, 60, 12},
        {105, 90, 75, 15}
    });
    
    auto result = matrix * vector;
    XCTAssertEqual(result, expected, "Product of two patrices OK");
}

- (void)testTransposeVectorToMatrix {
    const Vector<int> vector({0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
    const Matrix<int> expected({
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

@end
