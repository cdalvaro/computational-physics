//
//  VectorTests.m
//  Tests
//
//  Created by Carlos David on 11/06/2018.
//  Copyright © 2018 cdalvaro. All rights reserved.
//

#import <XCTest/XCTest.h>

#import "../../TestsTools.h"
#import "../../../computational-physics/math/containers/vector.hpp"

using namespace cda::math::containers;


@interface VectorTests : XCTestCase

@end

@implementation VectorTests

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
    [TestsTools setDefaultWorkingDirectory];
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

- (void)testComparisonBetweenTwoVectors {
    const Vector<double> vector1(10, 2);
    Vector<double> vector2(vector1);
    
    XCTAssert(vector1 == vector2, "Equality comparison OK");
    XCTAssert(!(vector1 != vector2), "Inequality comparison OK");
    
    const Vector<double> vector3(5, 2);
    XCTAssert(!(vector3 == vector1), "Equality comparison OK");
    XCTAssert(vector3 != vector1, "Inequality comparison OK");
    
    vector2.ones();
    XCTAssert(!(vector2 == vector1), "Equality comparison OK");
    XCTAssert(vector2 != vector1, "Inequality comparison OK");
}

- (void)testIsEmpty {
    const Vector<double> vector1(10, 2);
    XCTAssert(!vector1.IsEmpty(), "vector1 is not empty");
    
    const Vector<double> vector2;
    XCTAssert(vector2.IsEmpty(), "vector2 is empty");
    XCTAssertEqual(vector2.size(), 0, "The size of vector2 is 0");
}

- (void)testIsNull {
    const Vector<double> vector1(10, 0);
    XCTAssert(vector1.IsNull(), "vector1 is null");
    XCTAssert(!vector1.IsEmpty(), "vector1 is not empty");
    
    const Vector<double> vector2;
    XCTAssert(vector2.IsNull(), "vector2 is null");
    XCTAssert(vector2.IsEmpty(), "vector2 is empty");
    
    const Vector<double> vector3({0, 0, 0, 0, 1, 0, 0, 0, 0});
    XCTAssert(!vector3.IsNull(), "vector3 is not null");
    XCTAssert(!vector3.IsEmpty(), "vector3 is not empty");
}

- (void)testConstructors {
    
    // Default constructor
    const Vector<double> defaultVector;
    XCTAssert(defaultVector.IsEmpty(), "defaultVector is empty");
    XCTAssert(defaultVector.IsNull(), "defaultVector is null");
    
    // Constructor with size
    const Vector<double> vectorWithSize(10);
    XCTAssertEqual(vectorWithSize.size(), 10, "Constructor with size OK");
    XCTAssert(!vectorWithSize.IsEmpty() && !vectorWithSize.IsNull(), "vectorWithSize is not empty and is not null");
    
    // Constructor with size filling elements
    Vector<double> vectorWithElementsFilled(10, 5);
    XCTAssertEqual(vectorWithElementsFilled.size(), 10, "Constructor with size filling elements has size OK");
    XCTAssert(!vectorWithElementsFilled.IsEmpty() && !vectorWithElementsFilled.IsNull(), "vectorWithSize is not empty and is not null");
    for (size_t i = 0; i < vectorWithElementsFilled.size(); ++i) {
        XCTAssertEqual(vectorWithElementsFilled[i], 5, "Element has been filled OK");
    }
    
    // Copy constructor from another vector
    Vector<double> vectorFromCopy(vectorWithElementsFilled);
    XCTAssertEqual(vectorFromCopy, vectorWithElementsFilled, "Copy constructor OK");
    
    // Move constructor
    Vector<double> vectorFromMove(std::move(vectorWithElementsFilled));
    XCTAssertEqual(vectorFromMove, vectorFromCopy, "Move constructor OK");
    XCTAssert(vectorWithElementsFilled.IsEmpty(), "vectorWithElementsFilled is empty after had been moved");
    XCTAssert(vectorWithElementsFilled.IsNull(), "vectorWithElementsFilled is null after had been moved");
    
    // Constructor from array
    const double array[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    const Vector<double> vectorFromArray(array);
    XCTAssertEqual(vectorFromArray.size(), 10, "Constructor from array has right size");
    
    for (size_t i = 0; i < 10; ++i) {
        XCTAssertEqual(vectorFromArray[i], array[i], "Element has been copied OK");
    }
    
    const Vector<double> vectorFromArray2(array, array + 10);
    XCTAssertEqual(vectorFromArray2, vectorFromArray, "Vectro from array with size OK");
}

- (void)testCopy {
    const Vector<double> vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    
    Vector<double> vector_copy;
    vector_copy.copy(10, vector.begin());
    
    XCTAssertEqual(vector_copy, vector, "Copy of vectors OK");
}

- (void)testClear {
    Vector<double> vector(10, 1);
    
    vector.Clear();
    XCTAssert(vector.IsEmpty(), "Vector is empty");
    XCTAssert(vector.IsNull(), "Vector is null");
}

- (void)testIterators {
    Vector<double> vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    XCTAssertEqual(*vector.begin(), 1, "Begin iterator OK");
    XCTAssertEqual(*(vector.end() - 1), 10, "End iterator OK");
}

- (void)testAccessors {
    const Vector<double> const_vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    
    XCTAssertEqual(const_vector[2], 3, "Const accessor for reading OK");
    XCTAssertEqual(const_vector.at(5), 6, "Const at accessor for reading OK");
    XCTAssertThrows(const_vector.at(11), "Index out of bounds exception when accessing for reading OK");
    
    Vector<double> var_vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    XCTAssertEqual(var_vector[2], 3, "Const accessor for reading OK");
    XCTAssertEqual((var_vector[6] = 21), 21, "Var accesor after setting OK");
    XCTAssertEqual((var_vector.at(5) = 34), 34, "Var at accessor after setting OK");
    XCTAssertThrows(var_vector.at(11), "Index out of bounds exception when accessing for reading OK");
}

- (void)testGet {
    const Vector<double> vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    const Vector<double> expected1({5, 6, 7, 8, 9, 10});
    XCTAssertEqual(vector.get(4), expected1, "get without length OK");
    
    const Vector<double> expected2({2, 3, 4, 5, 6, 7, 8, 9});
    XCTAssertEqual(vector.get(1, 8), expected2, "get with length OK");
    
    XCTAssertThrows(vector.get(0, 11), "Index out of bounds");
    XCTAssertThrows(vector.get(3, 8), "Index out of bounds");
}

- (void)testSet {
    const Vector<double> original_vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    const Vector<double> mod_vector({2, 4, 6, 8});
    
    auto vector1 = original_vector;
    vector1.set(3, mod_vector);
    
    const Vector<double> expected1({1, 2, 3, 2, 4, 6, 8, 8, 9, 10});
    XCTAssertEqual(vector1, expected1, "Set without length OK");
    
    auto vector2 = original_vector;
    vector2.set(6, mod_vector, 3);
    
    const Vector<double> expected2({1, 2, 3, 4, 5, 6, 2, 4, 6, 10});
    XCTAssertEqual(vector2, expected2, "Set with length OK");
    
    XCTAssertThrows(vector2.set(8, vector2, 4), "Index out of bounds");
}

- (void)testAssign {
    const Vector<double> vector({ 0,  1,  2,  3 });
    
    Vector<double> test;
    test = vector;
    
    XCTAssertEqual(test, vector, "Assignment OK");
    
    test = test;
    XCTAssertEqual(test, vector, "Self assignment OK");
    
    Vector<double> test2;
    test2 = std::move(test);
    
    XCTAssertEqual(test2, vector, "Move assignment OK");
    XCTAssert(test.IsNull(), "test vector is null after had been moved");
    
    //Save the diagnostic state
#pragma clang diagnostic push
    
    // Ignore: Explicitly moving variable of type 'Matrix<double>' to itself
#pragma clang diagnostic ignored "-Wself-move"
    test2 = std::move(test2);
    XCTAssertEqual(test2, vector, "Move assignment does not destroy object if self matrix OK");
    
    //Restore the disgnostic state
#pragma clang diagnostic pop
}

- (void)testResize {
    Vector<double> vector({0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
    
    // Increase size
    vector.resize(15, true);
    XCTAssert(vector.size() == 15, "Vector size has been increased successfully");
    
    Boolean newElementsFilled = true;
    for (NSInteger i = 10; i < vector.size(); ++i) {
        if (vector[i] != 0) {
            newElementsFilled = false;
            break;
        }
    }
    
    XCTAssert(newElementsFilled, "New elements have been filled");
    
    // Resize to the same size
    auto previousVector(vector);
    vector.resize(vector.size());
    
    XCTAssert(previousVector == vector, "The vector has not changed");
    
    // Decrease size
    vector.resize(5);
    XCTAssert(vector.size() == 5, "Vector size has been decreased successfully");
    
    Boolean previousElementsEqual = true;
    for (NSInteger i = 0; i < vector.size(); ++i) {
        if (vector[i] != previousVector[i]) {
            previousElementsEqual = false;
            break;
        }
    }
    
    XCTAssert(previousElementsEqual, "Firsts elements still equal");
}

- (void)testFill {
    Vector<double> vector(10);
    const Vector<double> expected({2, 2, 2, 2, 2, 2, 2, 2, 2, 2});
    
    vector.fill(2);
    XCTAssertEqual(vector, expected, "Fill vector OK");
}

- (void)testOnes {
    const Vector<double> expected_ones(5, 1);
    XCTAssertEqual(Vector<double>::ones(5), expected_ones, "ones static method OK");
    
    Vector<double> result(5, 2);
    result.ones();
    XCTAssertEqual(result, expected_ones, "ones method OK");
}

- (void)testZero {
    const Vector<double> expected_zero(5, 0);
    XCTAssertEqual(Vector<double>::zero(5), expected_zero, "zero static method OK");
    
    Vector<double> result(5, 2);
    result.zero();
    XCTAssertEqual(result, expected_zero, "zero method OK");
}

- (void)testRandom {
    auto vector1 = Vector<double>::random(10, -3, 4);
    XCTAssert(!vector1.IsEmpty(), "Random vector is not empty");
    XCTAssert(!vector1.IsNull(), "Random vector is not null");
    XCTAssert(vector1.max_element() < 4, "The maximum element is smaller than 1");
    XCTAssert(vector1.min_element() >= -3, "The minimum element is greater or equal than 0");
    
    Vector<double> vector2(10);
    vector2.random();
    XCTAssert(!vector2.IsEmpty(), "Random vector is not empty");
    XCTAssert(!vector2.IsNull(), "Random vector is not null");
    XCTAssert(vector2.max_element() < 1, "The maximum element is smaller than 1");
    XCTAssert(vector2.min_element() >= 0, "The minimum element is greater or equal than 0");
}

- (void)testHasDuplicate {
    Vector<double> vector({1.5, 2.3, 3.6, 4.68, 5.6, 6.3, 7.34, 8.43, 9.43, 8.4363});
    XCTAssert(!vector.HasDuplicate(), "vector does not have duplicate elements");
    
    XCTAssert(vector.HasDuplicate(1E-02), "vector has duplicate with accuracy 1E-02");
    XCTAssert(!vector.HasDuplicate(1E-03), "vector does not have duplicate with accuracy 1E-03");
    
    vector[4] = 8.43;
    XCTAssert(vector.HasDuplicate(), "vector has duplicate");
}

- (void)testAdditionofVectors {
    Vector<double> vector1({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    const Vector<double> vector2({3, 5, 3, 5, 8, 3, 6, 7, 1, 5});
    const Vector<double> expected({4, 7, 6, 9, 13, 9, 13, 15, 10, 15});
    XCTAssertEqual(vector1 + vector2, expected, "The addition of two vectors is OK");
    
    vector1 += vector2;
    XCTAssertEqual(vector1, expected, "The addition of two vectors into the first one is OK");
    
    const Vector<double> vector3({1, 2, 3, 4, 5, 6});
    XCTAssertThrows(vector1 + vector3, "Vectors size must be equal");
    XCTAssertThrows(vector1 += vector3, "Vectors size must be equal");
}

- (void)testSubtractionOfVectors {
    Vector<double> vector1({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    const Vector<double> vector2({3, 5, 3, 5, 8, 3, 6, 7, 1, 5});
    const Vector<double> expected({-2, -3, 0, -1, -3, 3, 1, 1, 8, 5});
    XCTAssertEqual(vector1 - vector2, expected, "The subtraction of two vectors is OK");
    
    vector1 -= vector2;
    XCTAssertEqual(vector1, expected, "The subtraction of two vectors into the first one is OK");
    
    const Vector<double> vector3({1, 2, 3, 4, 5, 6});
    XCTAssertThrows(vector1 - vector3, "Vectors size must be equal");
    XCTAssertThrows(vector1 -= vector3, "Vectors size must be equal");
}

- (void)testProducts {
    Vector<double> vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    
    XCTAssertEqual(vector * vector, 385, "Product with another vector OK");
    
    const Vector<double> expected({2, 4, 6, 8, 10, 12, 14, 16, 18, 20});
    XCTAssertEqual(vector * 2.0, expected, "Product with scalar from the right OK");
    XCTAssertEqual(2.0 * vector, expected, "Product with scalar from the left OK");
    
    vector *= 2.0;
    XCTAssertEqual(vector, expected, "Product with scalar into itself OK");
    
    vector.resize(20);
    XCTAssertThrows(vector * expected, "Vectors must be of the same size");
    
    XCTAssertEqual(-vector, -1.0 * vector, "Negative vector OK");
}

- (void)testCrossProduct3D {
    const Vector<double> xAxis({1, 0, 0});
    const Vector<double> yAxis({0, 1, 0});
    const Vector<double> zAxis({0, 0, 1});
    
    XCTAssertEqual(xAxis.CrossProduct3D(yAxis), zAxis, "Cross product OK");
    
    const Vector<double> another_vector(4, 1);
    XCTAssertThrows(xAxis.CrossProduct3D(another_vector), "Vectors must be of the same size");
    XCTAssertThrows(another_vector.CrossProduct3D(another_vector), "Vector must be of size 3");
}

- (void)testDivisions {
    Vector<double> vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    const Vector<double> expected({0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5});
    
    XCTAssertEqual(vector / 2.0, expected, "Division with scalar OK");
    
    vector /= 2.0;
    XCTAssertEqual(vector, expected, "Division with scalar into itself OK");
}

- (void)testRemainder {
    Vector<double> vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    const Vector<double> expected({1, 2, 0, 1, 2, 0, 1, 2, 0, 1});
    XCTAssertEqual(vector % 3, expected, "Remainder OK");
    
    vector %= 3;
    XCTAssertEqual(vector, expected, "Remainder into itself OK");
}

- (void)testPowElements {
    const Vector<double> vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    XCTAssertEqual(vector.PowElements(0), Vector<double>::ones(10), "Power 0 OK");
    
    XCTAssertEqual(vector.PowElements(1), vector, "Power 1 OK");
    
    const Vector<double> expected({1, 8, 27, 64, 125, 216, 343, 512, 729, 1000});
    XCTAssertEqual(vector.PowElements(3), expected, "Power 3 OK");
}

- (void)testSqrt {
    const Vector<double> vector({1, 4, 9, 16, 25, 36, 49, 64, 81, 100});
    const Vector<double> expected({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    
    XCTAssertEqual(vector.Sqrt(), expected, "Sqrt method OK");
}

- (void)testSumAllElements {
    const Vector<double> vector1({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    XCTAssertEqual(vector1.sum(), 55, "SumAllElments OK");
}

- (void)testFindMethod {
    const Vector<double> vector({1, 2, 3, 4, 5, 6, 7, 4, 9, 10});
    
    auto it_value_4 = vector.Find(4);
    XCTAssertEqual(*it_value_4, 4, "Value 4 has been found");
    
    auto it_first_value_4 = vector.begin() + 3;
    XCTAssertEqual(it_value_4, it_first_value_4, "Found value 4 is the firstone");
    
    auto it_value_12 = vector.Find(12);
    XCTAssertEqual(it_value_12, vector.end(), "Value 12 is not present");
}

- (void)testMaximumAndMinimumElements {
    const Vector<double> vector({
         1.134,   0.001,   2.523,  -0.231,     0.321, -312353.123,
         5.213,   6.312,  -7.142,   8.243,     9.234,      21.426,
        10.123, -11.321,  12.213, -13.213,    14.231,       1.213,
        -0.024,   0.314,  17.143,  18.143, -1913.136,  523251.316,
         5.432,  -6.236,   7.342,   8.324,    -9.341,      21.341
    });
    
    XCTAssertEqual(vector.max_element(), 523251.316, "Maximum element OK");
    XCTAssertEqual(vector.abs_max_element(), 523251.316, "Absolute maximum element OK");
    XCTAssertEqual(vector.abs_max_element_with_sign(), 523251.316, "Absolute maximum element with sign OK");
    
    XCTAssertEqual(vector.min_element(), -312353.123, "Minimum element OK");
    XCTAssertEqual(vector.abs_min_element(), 0.001, "Absolute minimum element OK");
    XCTAssertEqual(vector.abs_min_element_with_sign(), 0.001, "Absolute minimum element with sign OK");
}

- (void)testSort {
    Vector<double> vector({2, 6, 4, 1, 7, 9, 10, 3, 5, 8});
    const Vector<double> expected({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    
    vector.sort();
    XCTAssertEqual(vector, expected, "sort method OK");
}

- (void)testNorms {
    const Vector<double> vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    
    // Square Norm
    auto result = vector.square_norm();
    double expected = 385;
    XCTAssertEqual(result, expected, "The square norm of the vector is OK");
    
    // Norm
    result = vector.norm();
    expected = std::sqrt(expected);
    XCTAssertEqual(result, expected, "The norm of the vector is OK");
}

- (void)testNormalizedVector {
    const Vector<double> vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    const auto expected = vector / std::sqrt(385);
    
    XCTAssertEqual(vector.normalized_vector(), expected, "Normalized vector OK");
}

- (void)testIfstreamOperator {
    
    std::ifstream file("data/math/containers/BigVector.csv", std::ios::in);
    
    // Big Matrix
    Vector<double> vector;
    file >> vector;
    file.close();
    
    Vector<double> expected_vector({
        306.956384885902, 292.352673102343, 287.916528879135, 285.99426058373, 284.989557035074,
        284.398942720766, 284.022457793086, 283.767793317283, 283.587526700195, 283.455258319253,
        283.35534128212, 283.278022793722, 283.216966971735, 283.167912596374, 283.127908123049,
        283.094856685579, 283.067235401401, 283.04391663485, 283.024050984977, 283.006988784374,
        282.992226275189, 282.979367978413, 282.968099918426, 282.958170265926, 282.949375139851,
        282.941548054634, 282.934551981049, 282.928273306183, 282.922617190522, 282.917503964648,
        282.912866307774, 282.908647020063, 282.904797250096, 282.901275074237, 282.898044350235,
        282.895073786212, 282.892336179983, 282.88980779401, 282.887467839031, 282.885298045301,
        282.883282304852, 282.881406371647, 282.879657609149, 282.878024776948, 282.876497849658,
        282.875067862659, 282.8737267802, 282.872467382255, 282.871283167131, 282.870168267386,
        282.869117377001, 282.868125688127, 282.867188835998, 282.866302850813, 282.865464115613,
        282.864669329309, 282.863915474141, 282.863199786996, 282.862519734047, 282.861872988298,
        282.861257409645, 282.860671027152, 282.860112023255, 282.859578719657, 282.859069564718,
        282.858583122154, 282.858118060898, 282.857673145985, 282.857247230341, 282.856839247386,
        282.856448204348, 282.856073176216, 282.855713300268, 282.855367771104, 282.855035836132,
        282.854716791472, 282.854409978216, 282.854114779028, 282.853830615032, 282.853556942976,
        282.853293252629, 282.853039064404, 282.852793927176, 282.852557416276, 282.852329131654,
        282.852108696186, 282.851895754111, 282.851689969602, 282.851491025432, 282.851298621759,
        282.851112474991, 282.850932316743, 282.850757892866, 282.850588962551, 282.850425297501,
        282.850266681154, 282.850112907969, 282.84996378276, 282.849819120076, 282.849678743623,
        282.849542485727, 282.849410186836, 282.849281695049, 282.849156865679, 282.84903556085,
        282.848917649113, 282.848803005091, 282.848691509145, 282.848583047059, 282.848477509753,
        282.848374793006, 282.848274797197, 282.848177427064, 282.848082591479, 282.847990203234,
        282.847900178838, 282.847812438331, 282.847726905107, 282.847643505743, 282.847562169844,
        282.847482829896, 282.847405421122, 282.847329881352, 282.847256150899, 282.847184172438,
        282.847113890897, 282.847045253352, 282.846978208925, 282.846912708693, 282.846848705594,
        282.846786154348, 282.846725011371, 282.846665234705, 282.846606783941, 282.846549620155,
        282.84649370584, 282.846439004845, 282.846385482316, 282.846333104645, 282.846281839411,
        282.846231655332, 282.84618252222, 282.846134410933, 282.84608729333, 282.846041142235,
        282.845995931392, 282.845951635433, 282.845908229838, 282.845865690904, 282.845823995712,
        282.845783122093, 282.845743048605, 282.845703754499, 282.845665219693, 282.845627424751,
        282.845590350852, 282.845553979771, 282.845518293855, 282.845483276002, 282.845448909639,
        282.845415178706, 282.845382067633, 282.845349561325, 282.845317645143, 282.845286304888,
        282.845255526787, 282.845225297475, 282.845195603982, 282.84516643372, 282.845137774466,
        282.845109614355, 282.845081941862, 282.845054745794, 282.845028015277, 282.845001739745,
        282.84497590893, 282.844950512853, 282.844925541813, 282.844900986377, 282.844876837374,
        282.844853085884, 282.844829723229, 282.844806740968, 282.844784130888, 282.844761884995,
        282.84473999551, 282.84471845486, 282.844697255672, 282.844676390766, 282.844655853151,
        282.844635636018, 282.844615732733, 282.844596136833, 282.844576842021, 282.844557842161,
        282.84453913127, 282.844520703519, 282.844502553224, 283.197374322374
    });
    
    XCTAssertEqual(vector, expected_vector, "Vector has been loaded properly.");
    
    std::ifstream unavailable_file("data/unavailable_file", std::ios::in);
    Matrix<double> matrix1;
    
    XCTAssertThrows(unavailable_file >> matrix1, "File does not exist OK");
    unavailable_file.close();
}

- (void)testOfstreamOperator {
    
    std::string expected_content;
    {
        std::ifstream expected_file;
        expected_file.open("data/math/containers/VectorTests-testOfstreamOperator.txt", std::ios::in);
        expected_content = std::string(std::istreambuf_iterator<char>(expected_file),
                                       std::istreambuf_iterator<char>());
        expected_file.close();
    }
    
    // --- Test
    std::stringstream test_output;
    std::streambuf *coutbuf = std::cout.rdbuf();
    
    // Redirect std::cout buffer
    std::cout.rdbuf(test_output.rdbuf());
    
    const Vector<double> vector1({21, 18, 15,  4, 42, 36, 30,  8, 84, 72, 63, 12});
    std::cout << vector1;
    
    // Restore std::cout buffer
    std::cout.rdbuf(coutbuf);
    test_output << vector1;
    
    XCTAssertEqual(test_output.str(), expected_content, "Ofstream operator OK");
}

@end
