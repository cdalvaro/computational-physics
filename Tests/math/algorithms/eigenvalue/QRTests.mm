//
//  QRTests.mm
//  Tests
//
//  Created by Carlos David on 15/08/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#import <XCTest/XCTest.h>

#import "../../../TestsTools.h"
#import "../../../../NumericalPDEs/math/algorithms/eigenvalue/qr.hpp"

using namespace cda::math::containers;
using namespace cda::math::algorithms::eigenvalue;


@interface QRTests : XCTestCase

@end

@implementation QRTests

const Matrix<double> test_matrix({
    {3,  2,  1,  2},
    {7,  8,  5,  1},
    {12, 11,  9,  8},
    {15, 14, 13, 12}
});

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
    [TestsTools setDefaultWorkingDirectory];
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

- (void)testQRDecomposition {
    QR<double> qr(test_matrix);
    
    const auto accuracy = 1E-4;
    
    // Matrix Q
    const Matrix<double> expected_q({
        {-0.1452,    0.5046,    0.6740,    0.5197},
        {-0.3388,   -0.8287,    0.4060,    0.1834},
        {-0.5807,    0.2128,    0.2813,   -0.7337},
        {-0.7259,    0.1156,   -0.5493,    0.3974}
    });
    
    auto result_q = qr.Q();
    
    for (auto it_result = result_q.Begin(), it_expected = expected_q.Begin();
         it_result != result_q.End(); ++it_result, ++it_expected) {
        XCTAssertEqualWithAccuracy(*it_result, *it_expected, accuracy, "Element of Q matrix is equal");
    }
    
    // Matrix R
    const Matrix<double> expected_r({
        {-20.6640,  -19.5509,   -16.5021,   -13.9857},
        {0.0000,     -1.6617,    -0.2213,     3.2698},
        {0.0000,      0.0000,    -1.9053,    -2.5873},
        {0.0000,      0.0000,     0.0000,     0.1223}
    });
    
    auto result_r = qr.R();
    
    for (auto it_result = result_r.Begin(), it_expected = expected_r.Begin();
         it_result != result_r.End(); ++it_result, ++it_expected) {
        XCTAssertEqualWithAccuracy(*it_result, *it_expected, accuracy, "Element of R matrix is equal");
    }
}

- (void)testQREigenValues {
    QR<double> qr(test_matrix);
    auto result_eigenvalues = qr.EigenValues();
    
    Vector<double> expected_eigenvalues({ 27.1593, 4.26358, 0.407643, 0.169475 });
    const auto accuracy = 1E-4;
    
    for (auto it_result = result_eigenvalues.Begin(), it_expected = expected_eigenvalues.Begin();
         it_result != result_eigenvalues.End(); ++it_result, ++it_expected) {
        XCTAssertEqualWithAccuracy(*it_result, *it_expected, accuracy, "Eigenvalue is OK");
    }
}

//- (void)testQREigenVectors {
//    QR<double> qr(test_matrix);
//    auto eigenvalues = qr.EigenValues();
//    
//    const Matrix<double> expected_eigenvectors({
//        {0.135474,   0.28509,   0.702763,   1},
//        { 0.17272,  -1.06963,   0.357511,   1},
//        {-2.76475,   3.10788,   -1.04857,   1},
//        {-2.00446,   2.46012,   -1.24656,   1}
//    });
//    
//    const auto accuracy = 1E-5;
//    
//    for (auto it = 0; it < eigenvalues.Size(); ++it) {
//        auto result = qr.EigenVector(eigenvalues.At(it));
//        auto expected = expected_eigenvectors.GetRowAsVector(it);
//        
//        for (auto it_result = result.Begin(), it_expected = expected.Begin();
//             it_result != result.End(); ++it_result, ++it_expected) {
//            XCTAssertEqualWithAccuracy(*it_result, *it_expected, accuracy, "Elements of eigen vector OK");
//        }
//    }
//}

- (void)testPerformanceExample {
    // This is an example of a performance test case.
    [self measureBlock:^{
        // Put the code you want to measure the time of here.
    }];
}

@end
