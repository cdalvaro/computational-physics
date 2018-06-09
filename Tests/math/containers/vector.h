//
//  vector.h
//  Wrapper for Vector class
//
//  Created by Carlos David on 09/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface VectorWrapper : NSObject

- (Boolean)comparisonBetweenTwoVectors;
- (Boolean)sumOfTwoVectors;
- (void)moveConstructorWithVectorSize:(NSInteger)size;

@end
