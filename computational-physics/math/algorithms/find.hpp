//
//  find.hpp
//  Computational Physics
//
//  Created by Carlos David on 23/08/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#pragma once

namespace cda {
    namespace math {
        namespace algorithms {
            namespace find {
                
                template <class InputIt, class T>
                constexpr InputIt Element(InputIt begin, InputIt const end, const T &value) {
                    for (; begin != end; ++begin) {
                        if (*begin == value) {
                            return begin;
                        }
                    }
                    return end;
                }
                
                template <typename T>
                T max_element(const T* const begin, const T* const end) {
                    auto max_element = *begin;
                    for (auto it = std::next(begin); it != end; ++it) {
                        if (*it > max_element) {
                            max_element = *it;
                        }
                    }
                    return max_element;
                }
                
                template <typename T>
                T abs_max_element(const T* const begin, const T* const end) {
                    T abs_it, max_element = std::abs(*begin);
                    for (auto it = std::next(begin); it != end; ++it) {
                        if ((abs_it = std::abs(*it)) > max_element) {
                            max_element = abs_it;
                        }
                    }
                    return max_element;
                }
                
                template <typename T>
                T AbsoluteMaximumElementWithSign(const T* const begin, const T* const end) {
                    auto max_element = *begin;
                    for (auto it = std::next(begin); it != end; ++it) {
                        if (std::abs(*it) > std::abs(max_element)) {
                            max_element = *it;
                        }
                    }
                    return max_element;
                }
                
                template <typename T>
                T MinimumElement(const T* const begin, const T* const end) {
                    auto min_element = *begin;
                    for (auto it = std::next(begin); it != end; ++it) {
                        if (*it < min_element) {
                            min_element = *it;
                        }
                    }
                    return min_element;
                }
                
                template <typename T>
                T AbsoluteMinimumElement(const T* const begin, const T* const end) {
                    T abs_it, min_element = std::abs(*begin);
                    for (auto it = std::next(begin); it != end; ++it) {
                        if ((abs_it = std::abs(*it)) < min_element) {
                            min_element = abs_it;
                        }
                    }
                    return min_element;
                }
                
                template <typename T>
                T AbsoluteMinimumElementWithSign(const T* const begin, const T* const end) {
                    auto min_element = *begin;
                    for (auto it = std::next(begin); it != end; ++it) {
                        if (std::abs(*it) < std::abs(min_element)) {
                            min_element = *it;
                        }
                    }
                    return min_element;
                }
                
            }
        }
    }
}
