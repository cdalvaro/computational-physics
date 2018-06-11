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
