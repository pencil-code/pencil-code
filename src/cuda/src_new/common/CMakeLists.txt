################################
##  CMakeLists.txt for utils  ##
################################

file (GLOB SRCS
     "*.cc")

add_library (common ${SRCS})
target_include_directories (common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})














































