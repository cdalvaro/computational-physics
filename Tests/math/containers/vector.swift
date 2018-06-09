//
//  vector.swift
//  Tests for Vector class
//
//  Created by Carlos David on 09/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

import XCTest

class vector: XCTestCase {
    
    var vectorWrapper = VectorWrapper()

    override func setUp() {
        super.setUp()
        // Put setup code here. This method is called before the invocation of each test method in the class.
    }

    override func tearDown() {
        super.tearDown()
        // Put teardown code here. This method is called after the invocation of each test method in the class.
    }
    
    func testComparisonBetweenTwoVectors(){
        XCTAssert(vectorWrapper.comparisonBetweenTwoVectors())
    }
    
    func testSumOfTwoVectors() {
        XCTAssert(vectorWrapper.sumOfTwoVectors())
    }

    func testPerformanceVectorMoveConstructor() {
        self.measure {
            for _ in 0...1_000 {
                vectorWrapper.moveConstructor(withVectorSize: 10_000)
            }
        }
    }

}
