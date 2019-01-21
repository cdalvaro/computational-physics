cmake_minimum_required( VERSION 3.5 )
project( ComputationalPhysics )

find_package( OpenGL REQUIRED )
find_package( GLUT REQUIRED )
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )

set( CMAKE_CXX_STANDARD 17 )

set( EXECUTABLE ComputationalPhysics )

include_directories( computational-physics )

file( GLOB_RECURSE SOURCE_FILES "computational-physics/*.cpp" )

add_executable( ${EXECUTABLE} ${SOURCE_FILES} )

target_link_libraries(
    ${EXECUTABLE}
    ${OPENGL_LIBRARIES}
    ${GLUT_LIBRARIES}
    pthread
)
