//
//  VectorTests.m
//  Tests
//
//  Created by Carlos David on 11/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#import <XCTest/XCTest.h>

#import "../../../NumericalPDEs/math/containers/vector.hpp"

using namespace cda::math::containers;

@interface VectorTests : XCTestCase

@end

@implementation VectorTests

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
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
    Vector<int> defaultVector;
    XCTAssert(defaultVector.Size() == 0, "Default constructor OK");
    
    // Constructor with size
    Vector<int> vectorWithSize(10);
    XCTAssert(vectorWithSize.Size() == 10, "Constructor with size OK");
    
    // Constructor with size filling elements
    Vector<int> vectorWithElementsFilled(10, 5);
    XCTAssert(vectorWithElementsFilled.Size() == 10, "Constructor with size filling elements has size OK");
    for (size_t i = 0; i < vectorWithElementsFilled.Size(); ++i) {
        XCTAssert(vectorWithElementsFilled[i] == 5, "Element has been filled OK");
    }
    
    // Copy constructor from another vector
    Vector<int> vectorFromCopy(vectorWithElementsFilled);
    XCTAssert(vectorFromCopy == vectorWithElementsFilled, "Copy constructor OK");
    
    // Move constructor
    Vector<int> vectorFromMove(std::move(vectorWithElementsFilled));
    XCTAssert(vectorFromMove == vectorFromCopy, "Move constructor OK");
    
    // Constructor from array
    const double array[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    Vector<double> vectorFromArray(array);
    XCTAssert(vectorFromArray.Size() == 10, "Constructor from array has right size");
    
    for (size_t i = 0; i < 10; ++i) {
        XCTAssert(vectorFromArray[i] == array[i], "Element has been copied OK");
    }
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
    
    XCTAssert(result == expected, "The sum of two vectors is OK");
}

- (void)testPerformanceVectorMoveConstructor {
    [self measureBlock:^{
        for (NSInteger i = 0; i < 1E+03; ++i) {
            Vector<double> vector(1E+04, 1);
            Vector<double> newVector(std::move(vector));
        }
    }];
}

@end
