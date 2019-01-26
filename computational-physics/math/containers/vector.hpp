//
//  vector.hpp
//  Class for vectorial operations
//  Vibration Normal Modes 2D
//
//  Created by Carlos David on 5/8/12.
//  Copyright (c) 2012 cdalvaro. All rights reserved.
//

#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "../algorithms/find.hpp"


namespace cda {
    namespace math {
        namespace containers {
    
            template <typename T>
            class Vector {
            private:
                size_t n;
                T *v, *it_end;
                
                void alloc_memory(const size_t &size) {
                    if (size == 0) {
                        std::free(v);
                        v = it_end = nullptr;
                    } else {
                        if (void *mem = std::realloc(v, size * sizeof(T))) {
                            v = static_cast<T *>(mem);
                            it_end = v + size;
                        } else {
                            throw std::bad_alloc();
                        }
                    }
                }
                
            public:
                
                typedef T value_type;
                
                /**
                 Class constructor
                 
                 @param size The size of the vector
                 */
                Vector(const size_t &size = 0) :
                n(size), v(nullptr), it_end(nullptr) {
                    alloc_memory(n);
                }
                
                /**
                 Class constructor
                 
                 @param size The size of the vector
                 */
                Vector(const size_t &size, const value_type &value) :
                n(size), v(nullptr), it_end(nullptr) {
                    alloc_memory(n);
                    fill(value);
                }
                
                /**
                 Copy constructor
                 
                 @param vector The source vector
                 */
                Vector(const Vector<value_type> &vector) :
                n(vector.n), v(nullptr), it_end(nullptr) {
                    alloc_memory(n);
                    std::copy(vector.begin(), vector.end(), this->begin());
                }
                
                /**
                 Move constructor
                 
                 @param vector The source vector
                 */
                Vector(Vector<value_type> &&vector) :
                n(vector.n), v(vector.v), it_end(vector.it_end) {
                    vector.n = 0;
                    vector.v = vector.it_end = nullptr;
                }
                
                template <size_t size>
                Vector(const value_type (& values)[size]) :
                n(size), v(nullptr), it_end(nullptr) {
                    alloc_memory(size);
                    std::copy(values, values + size, this->begin());
                }
                
                Vector(const value_type* const it_begin, const value_type* const it_end) :
                n(std::distance(it_begin, it_end)), v(nullptr), it_end(nullptr) {
                    alloc_memory(n);
                    std::copy(it_begin, it_end, this->begin());
                }
                
                /**
                 Class destructor
                 */
                ~Vector() {
                    if (v) {
                        std::free(v);
                        v = it_end = nullptr;
                        n = 0;
                    }
                }
                
                value_type *begin() const {
                    return v;
                }
                
                value_type *end() const {
                    return it_end;
                }
                
                /**
                 Change the size of the vector
                 
                 @param size The new size of the vector
                 @param fill If true and size is bigger than the previous size
                 the new extra elements will be set to 0
                 */
                void resize(const size_t &size, const bool &fill = true) {
                    if (size == n) {
                        return;
                    }
                    
                    alloc_memory(size);
                    if (fill && size > n) {
                        std::fill_n(v + n, size - n, static_cast<value_type>(0));
                    }
                    
                    n = size;
                }
                
                Vector<value_type> &operator=(const Vector<value_type> &vector) {
                    if (this != &vector) {
                        resize(vector.n, false);
                        std::copy(vector.begin(), vector.end(), this->begin());
                    }
                    return *this;
                }
                
                Vector<value_type> &operator=(Vector<value_type> &&vector) {
                    if (this != &vector) {
                        std::free(v);
                        v = vector.v;
                        it_end = vector.it_end;
                        n = vector.n;
                        
                        vector.v = vector.it_end = nullptr;
                        vector.n = 0;
                    }
                    return *this;
                }
                
                void copy(const size_t &size, const value_type* const array) {
                    resize(size);
                    std::copy(array, array + size, v);
                }
                
                bool operator==(const Vector<value_type> &vector) const {
                    if (this->n != vector.n) {
                        return false;
                    }
                    
                    auto it_vector = vector.begin();
                    for (auto it = this->begin(); it != this->end(); ++it, ++it_vector) {
                        if (*it != *it_vector) {
                            return false;
                        }
                    }
                    
                    return true;
                }
                
                bool operator!=(const Vector<value_type> &vector) const {
                    return !this->operator==(vector);
                }
                
                /**
                 Returns \p elements starting from element \p first_element
                 
                 @param first_element The first element to be recovered
                 @param elements The number of elements to be recovered
                 
                 @return A vector with size \p elements and the values starting from \p first_element
                 */
                Vector<value_type> get(const size_t &first_element, const size_t &elements) const {
                    if (n - first_element < elements) {
                        throw std::out_of_range("There are not enough elements inside the vector");
                    }
                    
                    auto it_begin = v + first_element;
                    auto it_end = it_begin + elements;
                    
                    return Vector<value_type>(it_begin, it_end);
                }
                
                /**
                 Returns a vector starting at position \p first_element
                 
                 @param first_element The position of the first element to be recovered
                 
                 @return A vector with all values starting from \p first_element
                 */
                Vector<value_type> get(const size_t &first_element) const {
                    return get(first_element, n - first_element);
                }
                
                /**
                 Copy \p elements from \p vector starting at position \p first_element
                 
                 @param first_element The position of the first element to be copied
                 @param vector The vector with the elements to be copied
                 @param elements The number of elements to copy
                 */
                void set(const size_t &first_element, const Vector<value_type> &vector, const size_t &elements) {
                    if (first_element + elements > n || elements > vector.n) {
                        throw std::out_of_range("Out of bounds");
                    }
                    std::copy(vector.v, vector.v + elements, v + first_element);
                }
                
                /**
                 Copy all elements of \p vector at position \p first_element
                 
                 @param first_element The position of the first element to be copied
                 @param vector The vector to copy
                 */
                void set(const size_t &first_element, const Vector<value_type> &vector) {
                    set(first_element, vector, vector.n);
                }
                
                size_t size() const {
                    return n;
                }
                
                value_type max_element() const {
                    return algorithms::find::max_element(begin(), end());
                }
                
                value_type abs_max_element() const {
                    return algorithms::find::abs_max_element(begin(), end());
                }
                
                value_type abs_max_element_with_sign() const {
                    return algorithms::find::abs_max_element_with_sign(begin(), end());
                }
                
                value_type min_element() const {
                    return algorithms::find::min_element(begin(), end());
                }
                
                value_type abs_min_element() const {
                    return algorithms::find::abs_min_element(begin(), end());
                }
                
                value_type abs_min_element_with_sign() const {
                    return algorithms::find::abs_min_element_with_sign(begin(), end());
                }
                
                value_type sum() const {
                    auto _sum = *begin();
                    for (auto it = begin() + 1; it != end(); ++it) {
                        _sum += *it;
                    }
                    return _sum;
                }
                
                const value_type &at(const size_t &element) const {
                    if (element >= n) {
                        throw std::out_of_range("Index out of bounds");
                    }
                    return v[element];
                }
                
                const value_type &operator[](const size_t &element) const {
                    return v[element];
                }
                
                value_type &at(const size_t &element) {
                    if (element >= n) {
                        throw std::out_of_range("Index out of bounds");
                    }
                    return v[element];
                }
                
                value_type &operator[](const size_t &element) {
                    return v[element];
                }
                
                double norm() const {
                    return std::sqrt(square_norm());
                }
                
                double square_norm() const {
                    return *this * *this;
                }
                
                /**
                 Returns a normalized vector
                 
                 @return The normalized vector
                 */
                Vector<value_type> normalized_vector() const {
                    return (* this) / norm();
                }
                
                void sort() {
                    std::sort(this->begin(), this->end(),
                              [](const value_type &value1, const value_type &value2) {
                                  return value1 < value2;
                              });
                }
                
                void fill(const value_type &value) {
                    std::fill(begin(), end(), value);
                }
                
                static Vector<value_type> zero(const size_t &size) {
                    return Vector<value_type>(size, 0);
                }
                
                void zero() {
                    fill(0);
                }
                
                static Vector<value_type> ones(const size_t &size) {
                    return Vector<value_type>(size, 1);
                }
                
                void ones() {
                    fill(1);
                }
                
                static Vector<value_type> random(const size_t &size, const value_type &min, const value_type &max) {
                    Vector<value_type> _random(size);
                    _random.random(min, max);
                    return _random;
                }
                
                void random(const value_type &min = 0, const value_type &max = 1) {
                    for (auto it = begin(); it != end(); ++it) {
                        *it = static_cast<value_type>(drand48() * max + min);
                    }
                }
                
                void clear() {
                    n = 0;
                    alloc_memory(n);
                }
                
                bool is_empty() const {
                    return n == 0;
                }
                
                bool is_null() const {
                    for (auto it = begin(); it != end(); ++it) {
                        if (*it != 0) {
                            return false;
                        }
                    }
                    return true;
                }
                
                bool has_duplicate(const value_type &accuracy) const {
                    value_type distance;
                    for (auto it1 = begin(); it1 != end(); ++it1) {
                        for (auto it2 = begin(); it2 != end(); ++it2) {
                            if (it1 != it2) {
                                distance = *it1 - *it2;
                                if (std::sqrt(distance * distance) < accuracy) {
                                    return true;
                                }
                            }
                        }
                    }
                    return false;
                }
                
                bool has_duplicate() const {
                    for (auto it1 = begin(); it1 != end(); ++it1) {
                        for (auto it2 = begin(); it2 != end(); ++it2) {
                            if (it1 != it2 && *it1 == *it2) {
                                return true;
                            }
                        }
                    }
                    return false;
                }
                
                value_type * find(const value_type &value) const {
                    return cda::math::algorithms::find::element(begin(), end(), value);
                }
                
                //  --- OPERATORS ---
                Vector<value_type> operator+(const Vector<value_type> &vector) const {
                    if (this->n != vector.n) {
                        throw std::logic_error("Unable to sum two vector of different size");
                    }
                    
                    Vector<value_type> new_vector(this->n);
                    
                    auto it_new = new_vector.begin();
                    auto it_vector = vector.begin();
                    for (auto it_this = this->begin(); it_this != this->end(); ++it_this, ++it_vector, ++it_new) {
                        *it_new = *it_this + *it_vector;
                    }
                    
                    return new_vector;
                }
                
                Vector<value_type> &operator+=(const Vector<value_type> &vector) {
                    if (this->n != vector.n) {
                        throw std::logic_error("Unable to sum two vector of different size");
                    }
                    
                    auto it_vector = vector.begin();
                    for (auto it_this = this->begin(); it_this != this->end(); ++it_this, ++it_vector) {
                        *it_this += *it_vector;
                    }
                    
                    return *this;
                }
                
                Vector<value_type> operator-(const Vector<value_type> &vector) const {
                    if (this->n != vector.n) {
                        throw std::logic_error("Unable to sum two vector of different sizes");
                    }
                    
                    Vector<value_type> new_vector(this->n);
                    
                    auto it_new = new_vector.begin();
                    auto it_vector = vector.begin();
                    for (auto it_this = this->begin(); it_this != this->end(); ++it_this, ++it_vector, ++it_new) {
                        *it_new = *it_this - *it_vector;
                    }
                    
                    return new_vector;
                }
                
                Vector<value_type> &operator-=(const Vector<value_type> &vector) {
                    if (this->n != vector.n) {
                        throw std::logic_error("Unable to sum two vector of different size");
                    }
                    
                    auto it_vector = vector.begin();
                    for (auto it_this = this->begin(); it_this != this->end(); ++it_this, ++it_vector) {
                        *it_this -= *it_vector;
                    }
                    
                    return *this;
                }
                
                Vector<value_type> operator*(const value_type &value) const {
                    Vector<value_type> new_vector(this->n);
                    auto it_new = new_vector.begin();
                    
                    for (auto it_this = this->begin(); it_this != this->end(); ++it_this, ++it_new) {
                        *it_new = *it_this * value;
                    }
                    
                    return new_vector;
                }
                
                Vector<value_type> &operator*=(const value_type &value) {
                    for (auto it = this->begin(); it != this->end(); ++it) {
                        *it *= value;
                    }
                    return *this;
                }
                
                Vector<value_type> operator/(const value_type &value) const {
                    Vector<value_type> new_vector(this->n);
                    auto it_new = new_vector.begin();
                    
                    for (auto it = this->begin(); it != this->end(); ++it, ++it_new) {
                        *it_new = *it / value;
                    }
                    
                    return new_vector;
                }
                
                Vector<value_type> &operator/=(const value_type &value) {
                    for (auto it = this->begin(); it != this->end(); ++it) {
                        *it /= value;
                    }
                    return *this;
                }
                
                Vector<value_type> pow(const size_t &power) const {
                    Vector<value_type> new_vector;
                    
                    switch (power) {
                        case 0:
                            new_vector = Vector<value_type>::ones(this->n);
                            break;
                        
                        default:
                            new_vector = *this;
                            value_type *it_new_vector;
                            for (size_t i = 2; i <= power; ++i) {
                                it_new_vector = new_vector.begin();
                                for (auto it_this = this->begin(); it_this != this->end(); ++it_this) {
                                    *it_new_vector++ *= *it_this;
                                }
                            }
                            break;
                    }
                    
                    return new_vector;
                }
                
                Vector<value_type> sqrt() const {
                    Vector<value_type> new_vector(this->n);
                    auto it_new_vector = new_vector.begin();
                    for (auto it_this = this->begin(); it_this != this->end(); ++it_this) {
                        *it_new_vector++ = std::sqrt(*it_this);
                    }
                    return new_vector;
                }
                
                template<typename Integer,
                         typename = std::enable_if<std::is_integral<Integer>::value>>
                Vector<value_type> operator%(const Integer &value) const {
                    Vector<value_type> new_vector(this->n);
                    auto it_new = new_vector.begin();
                    
                    for (auto it = this->begin(); it != this->end(); ++it, ++it_new) {
                        *it_new = static_cast<Integer>(*it) % value;
                    }
                    
                    return new_vector;
                }
                
                template<typename Integer,
                         typename = std::enable_if<std::is_integral<Integer>::value>>
                Vector<value_type> &operator%=(const Integer& value) {
                    for (auto it = this->begin(); it != this->end(); ++it) {
                        *it = static_cast<Integer>(*it) % value;
                    }
                    return *this;
                }
                
                Vector<value_type> operator-() const {
                    Vector<value_type> new_vector(this->n);
                    auto it_new = new_vector.begin();
                    
                    for (auto it = this->begin(); it != this->end(); ++it) {
                        *it_new++ = -(*it);
                    }
                    
                    return new_vector;
                }
                
                Vector<value_type> cross_product(const Vector<value_type> &vector) const {
                    if (this->n != vector.n || this->n != 3) {
                        throw std::logic_error("Both vectors must be of the same size, and size must be 3");
                    }
                    
                    Vector<value_type> new_vector(this->n);
                    new_vector[0] = this->operator[](1) * vector[2] - this->operator[](2) * vector[1];
                    new_vector[1] = this->operator[](2) * vector[0] - this->operator[](0) * vector[2];
                    new_vector[2] = this->operator[](0) * vector[1] - this->operator[](1) * vector[0];
                    
                    return new_vector;
                }
                
                value_type operator*(const Vector<value_type> &vector) const {
                    if (this->n != vector.n) {
                        throw std::logic_error("Both vectors must be of the same size");
                    }
                    
                    value_type value(0);
                    auto it_vector = vector.begin();
                    
                    for (auto it_this = this->begin(); it_this != this->end(); ++it_this, ++it_vector) {
                        value += *it_this * *it_vector;
                    }
                    
                    return value;
                }
                
            };
    
        } /* namespace math */
    } /* namespace containers */
} /* namespace containers */

//  --- MORE OPERATORS ---
template <typename ValueType, typename T2>
cda::math::containers::Vector<ValueType> operator*(const T2 &value, const cda::math::containers::Vector<ValueType> &vector) {
    
    cda::math::containers::Vector<ValueType> new_vector(vector.size());
    auto it_new = new_vector.begin();
    
    for (auto it = vector.begin(); it != vector.end(); ++it, ++it_new) {
        *it_new = *it * static_cast<ValueType>(value);
    }
    
    return new_vector;
}

template <typename ValueType>
void operator>>(std::istream &input,
                cda::math::containers::Vector<ValueType> &vector) {
    
    if (!input) {
        throw std::logic_error("Input is not avaiable");
    }
    
    size_t initial_size = 100;
    vector.resize(initial_size);
    
    char separator;
    std::string line, cell;
    size_t element = 0;
    if (std::getline(input, line)) {
        std::stringstream line_stream(line);
        while (line_stream >> vector[element]) {
            ++element;
            
            if (element >= vector.size()) {
                initial_size *= 2;
                vector.resize(vector.size() + initial_size);
            }
            
            line_stream >> separator;
        }
    }
    
    vector.resize(element);
}

template <typename ValueType>
std::ostream& operator<<(std::ostream &output,
                         const cda::math::containers::Vector<ValueType> &vector) {
    
    if (output.rdbuf() == std::cout.rdbuf()) {
        
        const size_t custom_width = 12;
        const size_t custom_precision = 5;
        
        output.width();
        output << std::fixed;
        output.fill(' ');
        output.precision(custom_precision);
        
        output << "[";
        output.width(custom_width);
        
        auto it_vector = vector.begin();
        output << *it_vector;
        for (it_vector = std::next(it_vector); it_vector < vector.end(); ++it_vector) {
            output << std::right << " ";
            output.width(custom_width);
            output << *it_vector;
        }
        output.width();
        output << "]" << std::endl;
        
        output.precision();
        
    } else {
        const char separator = ';';
        auto it_vector = vector.begin();
        output << *it_vector;
        for (it_vector = std::next(it_vector); it_vector < vector.end(); ++it_vector) {
            output << std::right << separator;
            output << *it_vector;
        }
        output << std::endl;
    }
    
    return output;
}
