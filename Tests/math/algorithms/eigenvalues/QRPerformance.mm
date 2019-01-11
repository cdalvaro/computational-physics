//
//  QRPerformance.mm
//  Tests
//
//  Created by Carlos David on 28/08/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#import <XCTest/XCTest.h>

#import "../../../TestsTools.h"
#import "../../../../computational-physics/math/algorithms/eigenvalues/qr.hpp"

using namespace cda::math::containers;
using namespace cda::math::algorithms::eigenvalues;


@interface QRPerformance : XCTestCase

@end

@implementation QRPerformance

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
    [TestsTools setDefaultWorkingDirectory];
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

- (void)testPerformanceEigenValues {
    
    std::ifstream file("data/math/algorithms/eigenvalues/Matrix-QRPerformance.csv", std::ios::in);
    
    Matrix<double> matrix;
    file >> matrix;
    file.close();
    
    [self measureBlock:^{
        QR<Matrix, double> qr(matrix, 1E-05, 1E+05);
        qr.EigenValues();
    }];
}

@end
