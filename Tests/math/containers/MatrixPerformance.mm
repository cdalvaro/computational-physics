//
//  MatrixPerformance.m
//  Tests
//
//  Created by Carlos David on 15/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#import <XCTest/XCTest.h>

#import "../../TestsTools.h"
#import "../../../NumericalPDEs/math/containers/matrix.hpp"

using namespace cda::math::containers;


@interface MatrixPerformance : XCTestCase

@end

@implementation MatrixPerformance

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
    [TestsTools setDefaultWorkingDirectory];
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

- (void)testPerformanceMatrixMoveConstructor {
    [self measureBlock:^{
        Matrix<double> matrix(1E+04, 1E+04, 1);
        Matrix<double> newMatrix(std::move(matrix));
    }];
}

- (void)testPerformanceMatrixProduct {
    [self measureBlock:^{
        Matrix<double> matrix(200, 200, 1);
        auto new_matrix = matrix * matrix;
    }];
}

- (void)testPerformanceLoadMatrixFromFile {
    [self measureBlock:^{
        std::ifstream file("data/math/containers/BigMatrix.csv", std::ios::in);
        Matrix<double> matrix;
        file >> matrix;
        file.close();
    }];
}

@end
