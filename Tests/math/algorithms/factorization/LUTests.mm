//
//  LUTests.m
//  Tests
//
//  Created by Carlos David on 16/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#import <XCTest/XCTest.h>

#import "../../../TestsTools.h"
#import "../../../../NumericalPDEs/math/algorithms/factorization/lu.hpp"

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

- (void)testLUMethods {
    Matrix<double> matrix_test({
        {3,  2,  1,  2},
        {7,  6,  5,  1},
        {12, 10,  9,  8},
        {15, 14, 13, 12}
    });
    
    LU<double> lu(matrix_test);
    
    const auto accuracy = 1E-13;
    
    // Matrix U
    Matrix<double> expected_u({
        {3.0,     2.0,     1.0,       2.0},
        {0.0, 4.0/3.0, 8.0/3.0, -11.0/3.0},
        {0.0,     0.0,     1.0,  11.0/2.0},
        {0.0,     0.0,     0.0,      13.0}
    });
    
    auto result_u = lu.U();
    
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
    
    auto result_l = lu.L();
    
    for (auto it_result = result_l.Begin(), it_expected = expected_l.Begin();
         it_result != result_l.End(); ++it_result, ++it_expected) {
        XCTAssertEqualWithAccuracy(*it_result, *it_expected, accuracy, "Element of matrix L is equal");
    }
    
}

@end
