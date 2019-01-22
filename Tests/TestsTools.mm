//
//  TestsTools.m
//  Tests
//
//  Created by Carlos David on 15/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "TestsTools.h"


@implementation TestsTools

+ (void)setDefaultWorkingDirectory {
    
    NSURL *testsPath = [[NSURL fileURLWithPath:@__FILE__] URLByDeletingLastPathComponent];
    if (![[NSFileManager defaultManager] changeCurrentDirectoryPath:[testsPath path]]) {
        @throw [NSException exceptionWithName:@"UnableToChangeWorkingDirectory" reason:NULL userInfo:NULL];
    }
    
}

+ (BOOL)compareMatrix: (cda::math::containers::Matrix<double>) compare
         withExpected: (cda::math::containers::Matrix<double>) expected
         whitAccuracy: (double) accuracy {
    
    if (compare.Rows() != expected.Rows() || compare.Columns() != expected.Columns()) {
        return NO;
    }
    
    for (auto it_compare = compare.begin(), it_expected = expected.begin();
         it_compare != compare.end(); ++it_compare, ++it_expected) {
        
        const double distance = *it_compare - *it_expected;
        if (std::sqrt(distance * distance) >= accuracy) {
            return NO;
        }
    }
    
    return YES;
}

+ (BOOL)compareVector: (cda::math::containers::Vector<double>) compare
         withExpected: (cda::math::containers::Vector<double>) expected
         whitAccuracy: (double) accuracy {
    
    if (compare.size() != expected.size()) {
        return NO;
    }
    
    for (auto it_compare = compare.begin(), it_expected = expected.begin();
         it_compare != compare.end(); ++it_compare, ++it_expected) {
        
        const double distance = *it_compare - *it_expected;
        if (std::sqrt(distance * distance) >= accuracy) {
            return NO;
        }
    }
    
    return YES;
}

@end
