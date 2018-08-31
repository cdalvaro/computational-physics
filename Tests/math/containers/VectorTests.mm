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
    Vector<double> vector1(10, 2);
    Vector<double> vector2(vector1);
    
    XCTAssert(vector1 == vector2, "Equality comparison OK");
    XCTAssert(!(vector1 != vector2), "Inequality comparison OK");
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
}

- (void)testIterators {
    Vector<double> vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    XCTAssertEqual(*vector.Begin(), 1, "Begin iterator OK");
    XCTAssertEqual(*(vector.End() - 1), 10, "End iterator OK");
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
    
    // Resize to the main size
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

- (void)testSumOfTwoVectors {
    Vector<double> vector1(10);
    for (size_t i = 0; i < vector1.Size(); ++i) {
        vector1[i] = i;
    }
    
    Vector<double> vector2(10);
    for (size_t i = 0; i < vector2.Size(); ++i) {
        vector2[i] = 10 + i;
    }
    
    Vector<double> expected(10);
    for (size_t i = 0; i < expected.Size(); ++i) {
        expected[i] = vector1[i] + vector2[i];
    }
    
    auto result = vector1 + vector2;
    
    XCTAssertEqual(result, expected, "The sum of two vectors is OK");
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
