//
//  QRTests.mm
//  Tests
//
//  Created by Carlos David on 15/08/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#import <XCTest/XCTest.h>

#import "../../../TestsTools.h"
#import "../../../../computational-physics/math/algorithms/eigenvalues/qr.hpp"

using namespace cda::math::containers;
using namespace cda::math::algorithms::eigenvalues;


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

- (void)testConstructor {
    Matrix<double> matrix_test({
        {3,  2,  1,  2},
        {7,  6,  5,  1},
        {12, 10,  9,  8},
        {15, 14, 13, 12}
    });
    
    XCTAssertNoThrow(QR<Matrix>(matrix_test), "QR matrix constructor OK");
    XCTAssertThrows(QR<Matrix>(Matrix<double>(4, 3)), "QR matrix is not valid for non-square matrices");
}

- (void)testQRDecomposition {
    QR<Matrix> qr(test_matrix);
    
    const auto accuracy = 1E-4;
    
    // Q Matrix
    const Matrix<double> expected_q({
        {-0.1452,    0.5046,    0.6740,    0.5197},
        {-0.3388,   -0.8287,    0.4060,    0.1834},
        {-0.5807,    0.2128,    0.2813,   -0.7337},
        {-0.7259,    0.1156,   -0.5493,    0.3974}
    });
    
    XCTAssert([TestsTools compareMatrix:qr.q()
                           withExpected:expected_q
                           whitAccuracy:accuracy],
              "Q matrix OK");
    
    // R Matrix
    const Matrix<double> expected_r({
        {-20.6640,  -19.5509,   -16.5021,   -13.9857},
        {0.0000,     -1.6617,    -0.2213,     3.2698},
        {0.0000,      0.0000,    -1.9053,    -2.5873},
        {0.0000,      0.0000,     0.0000,     0.1223}
    });
    
    XCTAssert([TestsTools compareMatrix:qr.r()
                           withExpected:expected_r
                           whitAccuracy:accuracy],
              "R matrix OK");
}

- (void)testQREigenValues {
    QR<Matrix> qr(test_matrix);
    
    const Vector<double> expected_eigenvalues({ 27.1593, 4.26358, 0.407643, 0.169475 });
    const auto accuracy = 1E-4;
    
    XCTAssert([TestsTools compareVector:qr.eigen_values()
                           withExpected:expected_eigenvalues
                           whitAccuracy:accuracy],
              "Eigenvalues OK");
    
    XCTAssert([TestsTools compareVector:qr.eigen_values()
                           withExpected:expected_eigenvalues
                           whitAccuracy:accuracy],
              "Eigenvalues cache OK");
}

- (void)testQREigenVectorsOneByOne {
    QR<Matrix> qr(test_matrix);
    auto eigenvalues = qr.eigen_values();
    
    const Matrix<double> expected_eigenvectors({
        { 0.13547,   0.28509,    0.70276,   1},
        { 0.17272,  -1.06963,    0.35751,   1},
        {-2.76475,   3.10788,   -1.04857,   1},
        {-2.00446,   2.46012,   -1.24656,   1}
    });
    
    const auto accuracy = 1E-5;
    
    for (auto it = 0; it < eigenvalues.size(); ++it) {
        XCTAssert([TestsTools compareVector:qr.EigenVector(eigenvalues.at(it))
                               withExpected:expected_eigenvectors.get_row_as_vector(it)
                               whitAccuracy:accuracy],
                  @"Eigenvector OK for eigenvalue %f", eigenvalues.at(it));
        
        XCTAssert([TestsTools compareVector:qr.EigenVector(eigenvalues.at(it))
                               withExpected:expected_eigenvectors.get_row_as_vector(it)
                               whitAccuracy:accuracy],
                  @"Eigenvector cache OK for eigenvalue %f", eigenvalues.at(it));
    }
}

- (void)testQREigenVectorsAtOnce {
    QR<Matrix> qr(test_matrix);
    auto eigenvalues = qr.eigen_values();
    auto eigenvectors = qr.EigenVectors();
    
    const Matrix<double> expected_eigenvectors({
        { 0.13547,   0.28509,    0.70276,   1},
        { 0.17272,  -1.06963,    0.35751,   1},
        {-2.76475,   3.10788,   -1.04857,   1},
        {-2.00446,   2.46012,   -1.24656,   1}
    });
    
    const auto accuracy = 1E-5;
    
    for (auto it = 0; it < eigenvalues.size(); ++it) {
        XCTAssert([TestsTools compareVector:eigenvectors.at(eigenvalues.at(it))
                               withExpected:expected_eigenvectors.get_row_as_vector(it)
                               whitAccuracy:accuracy],
                  @"Eigenvector OK for eigenvalue %f", eigenvalues.at(it));
    }
}

@end
