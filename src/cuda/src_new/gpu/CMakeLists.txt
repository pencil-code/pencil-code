################################
##  CMakeLists.txt for utils  ##
################################

file (GLOB SRCS
     "*.cc")

add_library (gpu ${SRCS})
target_include_directories (gpu PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})














































