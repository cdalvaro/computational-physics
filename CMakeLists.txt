cmake_minimum_required( VERSION 3.9.0 )
project( NumericalPDEs )

find_package( OpenGL )
find_package( GLUT )

set( CMAKE_CXX_STANDARD 17 )

set( EXECUTABLE GraphicalPDEs )

include_directories( NumericalPDEs )

file( GLOB_RECURSE SOURCE_FILES "NumericalPDEs/*.cpp" )

add_executable( ${EXECUTABLE} ${SOURCE_FILES} )

target_link_libraries(
    ${EXECUTABLE}
    ${OPENGL_gl_LIBRARY}
    ${GLUT_LIBRARIES} )

