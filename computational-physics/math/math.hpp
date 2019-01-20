//
//  math.hpp
//  Computational Physics
//
//  Created by Carlos David on 17/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#pragma once

#include <type_traits>


namespace cda {
    namespace math {
        
        template <typename T> inline constexpr
        T signum(const T &x, std::false_type is_signed) {
            return static_cast<T>(0) < x;
        }
        
        template <typename T> inline constexpr
        T signum(const T &x, std::true_type is_signed) {
            return (static_cast<T>(0) < x) - (x < static_cast<T>(0));
        }
        
        template <typename T> inline constexpr
        T signum(const T &x) {
            return signum(x, std::is_signed<T>());
        }
        
    }
}
