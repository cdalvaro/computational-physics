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

@end
