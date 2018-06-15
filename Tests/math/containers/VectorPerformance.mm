//
//  VectorPerformance.m
//  Tests
//
//  Created by Carlos David on 15/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#import <XCTest/XCTest.h>

#import "../../TestsTools.h"
#import "../../../NumericalPDEs/math/containers/vector.hpp"

using namespace cda::math::containers;


@interface VectorPerformance : XCTestCase

@end

@implementation VectorPerformance

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
    [TestsTools setDefaultWorkingDirectory];
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
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
