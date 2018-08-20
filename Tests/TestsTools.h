//
//  TestsTools.h
//  NumericalPDEs
//
//  Created by Carlos David on 15/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#ifndef TestsTools_h
#define TestsTools_h

#import "../NumericalPDEs/math/containers.hpp"


@interface TestsTools : NSObject

+ (void)setDefaultWorkingDirectory;

+ (BOOL)compareMatrix: (cda::math::containers::Matrix<double>) compare
         withExpected: (cda::math::containers::Matrix<double>) expected
         whitAccuracy: (double) accuracy;

@end

#endif /* Tools_h */
