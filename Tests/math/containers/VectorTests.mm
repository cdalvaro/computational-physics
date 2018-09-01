//
//  VectorTests.m
//  Tests
//
//  Created by Carlos David on 11/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#import <XCTest/XCTest.h>

#import "../../TestsTools.h"
#import "../../../NumericalPDEs/math/containers/vector.hpp"

using namespace cda::math::containers;


@interface VectorTests : XCTestCase

@end

@implementation VectorTests

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
    [TestsTools setDefaultWorkingDirectory];
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

- (void)testComparisonBetweenTwoVectors {
    const Vector<double> vector1(10, 2);
    Vector<double> vector2(vector1);
    
    XCTAssert(vector1 == vector2, "Equality comparison OK");
    XCTAssert(!(vector1 != vector2), "Inequality comparison OK");
    
    const Vector<double> vector3(5, 2);
    XCTAssert(!(vector3 == vector1), "Equality comparison OK");
    XCTAssert(vector3 != vector1, "Inequality comparison OK");
    
    vector2.Ones();
    XCTAssert(!(vector2 == vector1), "Equality comparison OK");
    XCTAssert(vector2 != vector1, "Inequality comparison OK");
}

- (void)testConstructors {
    
    // Default constructor
    const Vector<double> defaultVector;
    XCTAssertEqual(defaultVector.Size(), 0, "Default constructor OK");
    XCTAssert(defaultVector.IsEmpty(), "defaultVector is empty");
    XCTAssert(defaultVector.IsNull(), "defaultVector is null");
    
    // Constructor with size
    const Vector<double> vectorWithSize(10);
    XCTAssertEqual(vectorWithSize.Size(), 10, "Constructor with size OK");
    XCTAssert(!vectorWithSize.IsEmpty() && !vectorWithSize.IsNull(), "vectorWithSize is not empty and is not null");
    
    // Constructor with size filling elements
    Vector<double> vectorWithElementsFilled(10, 5);
    XCTAssertEqual(vectorWithElementsFilled.Size(), 10, "Constructor with size filling elements has size OK");
    XCTAssert(!vectorWithElementsFilled.IsEmpty() && !vectorWithElementsFilled.IsNull(), "vectorWithSize is not empty and is not null");
    for (size_t i = 0; i < vectorWithElementsFilled.Size(); ++i) {
        XCTAssertEqual(vectorWithElementsFilled[i], 5, "Element has been filled OK");
    }
    
    // Copy constructor from another vector
    Vector<double> vectorFromCopy(vectorWithElementsFilled);
    XCTAssertEqual(vectorFromCopy, vectorWithElementsFilled, "Copy constructor OK");
    
    // Move constructor
    Vector<double> vectorFromMove(std::move(vectorWithElementsFilled));
    XCTAssertEqual(vectorFromMove, vectorFromCopy, "Move constructor OK");
    XCTAssert(vectorWithElementsFilled.IsEmpty(), "vectorWithElementsFilled is empty after had been moved");
    XCTAssert(vectorWithElementsFilled.IsNull(), "vectorWithElementsFilled is null after had been moved");
    
    // Constructor from array
    const double array[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    const Vector<double> vectorFromArray(array);
    XCTAssertEqual(vectorFromArray.Size(), 10, "Constructor from array has right size");
    
    for (size_t i = 0; i < 10; ++i) {
        XCTAssertEqual(vectorFromArray[i], array[i], "Element has been copied OK");
    }
    
    const Vector<double> vectorFromArray2(array, array + 10);
    XCTAssertEqual(vectorFromArray2, vectorFromArray, "Vectro from array with size OK");
}

- (void)testCopy {
    const Vector<double> vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    
    Vector<double> vector_copy;
    vector_copy.Copy(10, vector.Begin());
    
    XCTAssertEqual(vector_copy, vector, "Copy of vectors OK");
}

- (void)testIterators {
    Vector<double> vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    XCTAssertEqual(*vector.Begin(), 1, "Begin iterator OK");
    XCTAssertEqual(*(vector.End() - 1), 10, "End iterator OK");
}

- (void)testAccessors {
    const Vector<double> const_vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    
    XCTAssertEqual(const_vector[2], 3, "Const accessor for reading OK");
    XCTAssertEqual(const_vector.At(5), 6, "Const At accessor for reading OK");
    XCTAssertThrows(const_vector.At(11), "Index out of bounds exception when accessing for reading OK");
    
    Vector<double> var_vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    XCTAssertEqual(var_vector[2], 3, "Const accessor for reading OK");
    XCTAssertEqual((var_vector[6] = 21), 21, "Var accesor after setting OK");
    XCTAssertEqual((var_vector.At(5) = 34), 34, "Var At accessor after setting OK");
    XCTAssertThrows(var_vector.At(11), "Index out of bounds exception when accessing for reading OK");
}

- (void)testGet {
    const Vector<double> vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    const Vector<double> expected1({5, 6, 7, 8, 9, 10});
    XCTAssertEqual(vector.Get(4), expected1, "Get without length OK");
    
    const Vector<double> expected2({2, 3, 4, 5, 6, 7, 8, 9});
    XCTAssertEqual(vector.Get(1, 8), expected2, "Get with length OK");
    
    XCTAssertThrows(vector.Get(0, 11), "Index out of bounds");
    XCTAssertThrows(vector.Get(3, 8), "Index out of bounds");
}

- (void)testSet {
    const Vector<double> original_vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    const Vector<double> mod_vector({2, 4, 6, 8});
    
    auto vector1 = original_vector;
    vector1.Set(3, mod_vector);
    
    const Vector<double> expected1({1, 2, 3, 2, 4, 6, 8, 8, 9, 10});
    XCTAssertEqual(vector1, expected1, "Set without length OK");
    
    auto vector2 = original_vector;
    vector2.Set(6, mod_vector, 3);
    
    const Vector<double> expected2({1, 2, 3, 4, 5, 6, 2, 4, 6, 10});
    XCTAssertEqual(vector2, expected2, "Set with length OK");
    
    XCTAssertThrows(vector2.Set(8, vector2, 4), "Index out of bounds");
}

- (void)testAssign {
    const Vector<double> vector({ 0,  1,  2,  3 });
    
    Vector<double> test;
    test = vector;
    
    XCTAssertEqual(test, vector, "Assignment OK");
    
    test = test;
    XCTAssertEqual(test, vector, "Self assignment OK");
    
    Vector<double> test2;
    test2 = std::move(test);
    
    XCTAssertEqual(test2, vector, "Move assignment OK");
    XCTAssert(test.IsNull(), "test vector is null after had been moved");
    
    //Save the diagnostic state
#pragma clang diagnostic push
    
    // Ignore: Explicitly moving variable of type 'Matrix<double>' to itself
#pragma clang diagnostic ignored "-Wself-move"
    test2 = std::move(test2);
    XCTAssertEqual(test2, vector, "Move assignment does not destroy object if self matrix OK");
    
    //Restore the disgnostic state
#pragma clang diagnostic pop
}

- (void)testResize {
    Vector<double> vector({0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
    
    // Increase size
    vector.Resize(15, true);
    XCTAssert(vector.Size() == 15, "Vector size has been increased successfully");
    
    Boolean newElementsFilled = true;
    for (NSInteger i = 10; i < vector.Size(); ++i) {
        if (vector[i] != 0) {
            newElementsFilled = false;
            break;
        }
    }
    
    XCTAssert(newElementsFilled, "New elements have been filled");
    
    // Resize to the same size
    auto previousVector(vector);
    vector.Resize(vector.Size());
    
    XCTAssert(previousVector == vector, "The vector has not changed");
    
    // Decrease size
    vector.Resize(5);
    XCTAssert(vector.Size() == 5, "Vector size has been decreased successfully");
    
    Boolean previousElementsEqual = true;
    for (NSInteger i = 0; i < vector.Size(); ++i) {
        if (vector[i] != previousVector[i]) {
            previousElementsEqual = false;
            break;
        }
    }
    
    XCTAssert(previousElementsEqual, "Firsts elements still equal");
}

- (void)testOnes {
    const Vector<double> expected_ones(5, 1);
    XCTAssertEqual(Vector<double>::Ones(5), expected_ones, "Ones static method OK");
    
    Vector<double> result(5, 2);
    result.Ones();
    XCTAssertEqual(result, expected_ones, "Ones method OK");
}

- (void)testZero {
    const Vector<double> expected_zero(5, 0);
    XCTAssertEqual(Vector<double>::Zero(5), expected_zero, "Zero static method OK");
    
    Vector<double> result(5, 2);
    result.Zero();
    XCTAssertEqual(result, expected_zero, "Zero method OK");
}

- (void)testAdditionofVectors {
    Vector<double> vector1({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    const Vector<double> vector2({3, 5, 3, 5, 8, 3, 6, 7, 1, 5});
    const Vector<double> expected({4, 7, 6, 9, 13, 9, 13, 15, 10, 15});
    XCTAssertEqual(vector1 + vector2, expected, "The addition of two vectors is OK");
    
    vector1 += vector2;
    XCTAssertEqual(vector1, expected, "The addition of two vectors into the first one is OK");
    
    const Vector<double> vector3({1, 2, 3, 4, 5, 6});
    XCTAssertThrows(vector1 + vector3, "Vectors size must be equal");
    XCTAssertThrows(vector1 += vector3, "Vectors size must be equal");
}

- (void)testSubtractionOfVectors {
    Vector<double> vector1({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    const Vector<double> vector2({3, 5, 3, 5, 8, 3, 6, 7, 1, 5});
    const Vector<double> expected({-2, -3, 0, -1, -3, 3, 1, 1, 8, 5});
    XCTAssertEqual(vector1 - vector2, expected, "The subtraction of two vectors is OK");
    
    vector1 -= vector2;
    XCTAssertEqual(vector1, expected, "The subtraction of two vectors into the first one is OK");
    
    const Vector<double> vector3({1, 2, 3, 4, 5, 6});
    XCTAssertThrows(vector1 - vector3, "Vectors size must be equal");
    XCTAssertThrows(vector1 -= vector3, "Vectors size must be equal");
}

- (void)testProducts {
    Vector<double> vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    
    XCTAssertEqual(vector * vector, 385, "Product with another vector OK");
    
    const Vector<double> expected({2, 4, 6, 8, 10, 12, 14, 16, 18, 20});
    XCTAssertEqual(vector * 2.0, expected, "Product with scalar from the right OK");
    XCTAssertEqual(2.0 * vector, expected, "Product with scalar from the left OK");
    
    vector *= 2.0;
    XCTAssertEqual(vector, expected, "Product with scalar into itself OK");
    
    vector.Resize(20);
    XCTAssertThrows(vector * expected, "Vectors must be of the same size");
}

- (void)testSumAllElements {
    const Vector<double> vector1({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    XCTAssertEqual(vector1.SumAllEments(), 55, "SumAllElments OK");
}

- (void)testMaximumAndMinimumElements {
    const Vector<double> vector({
         1.134,   0.001,   2.523,  -0.231,     0.321, -312353.123,
         5.213,   6.312,  -7.142,   8.243,     9.234,      21.426,
        10.123, -11.321,  12.213, -13.213,    14.231,       1.213,
        -0.024,   0.314,  17.143,  18.143, -1913.136,  523251.316,
         5.432,  -6.236,   7.342,   8.324,    -9.341,      21.341
    });
    
    XCTAssertEqual(vector.MaximumElement(), 523251.316, "MaximumElement OK");
    XCTAssertEqual(vector.AbsoluteMaximumElement(), 523251.316, "AbsoluteMaximumElement OK");
    XCTAssertEqual(vector.AbsoluteMaximumElementWithSign(), 523251.316, "AbsoluteMaximumElementWithSign OK");
    
    XCTAssertEqual(vector.MinimumElement(), -312353.123, "MinimumElement OK");
    XCTAssertEqual(vector.AbsoluteMinimumElement(), 0.001, "AbsoluteMinimumElement OK");
    XCTAssertEqual(vector.AbsoluteMinimumElementWithSign(), 0.001, "AbsoluteMinimumElementWithSign OK");
}

- (void)testNorms {
    const Vector<double> vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    
    // Square Norm
    auto result = vector.SquareNorm();
    double expected = 385;
    XCTAssertEqual(result, expected, "The square norm of the vector is OK");
    
    // Norm
    result = vector.Norm();
    expected = std::sqrt(expected);
    XCTAssertEqual(result, expected, "The norm of the vector is OK");
}

@end
