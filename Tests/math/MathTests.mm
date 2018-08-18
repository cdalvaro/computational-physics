//
//  MathTests.mm
//  Tests
//
//  Created by Carlos David on 18/08/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#import <XCTest/XCTest.h>

#import "../TestsTools.h"
#import "../../NumericalPDEs/math/math.hpp"

using namespace cda::math;


@interface MathTests : XCTestCase

@end

@implementation MathTests

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
    [TestsTools setDefaultWorkingDirectory];
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

- (void)testSignum {
    XCTAssertEqual(signum(-3), -1, "Number is negative");
    XCTAssertEqual(signum(4), 1, "Number is positive");
    XCTAssertEqual(signum(0), 0, "Zero is 0");
    
    XCTAssertEqual(signum(static_cast<unsigned int>(4)), 1, "Unsigned number is positive");
}

@end
