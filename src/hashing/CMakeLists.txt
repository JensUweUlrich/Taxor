cmake_minimum_required(VERSION 3.16)

file(GLOB HASH_SRC "*.cpp" "*.c")
file(GLOB HASH_HDRS "*.hpp" "*.h")


configure_file(CMakeLists.txt.in ${CMAKE_CURRENT_BINARY_DIR}/hashing-download/CMakeLists.txt)

execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hashing-download)
if(result)
  message(FATAL_ERROR "CMake step for ankerl failed: ${result}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build . --config Release
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hashing-download )


add_library (hashing STATIC ${HASH_SRC} ${HASH_HDRS})

set(unordered_dense_DIR "${CMAKE_CURRENT_BINARY_DIR}/ankerl/lib${LIBDIRSUFFIX}/cmake/unordered_dense")
find_package(unordered_dense CONFIG REQUIRED)


# include SafeQueue and StopClock
include_directories("${CMAKE_CURRENT_BINARY_DIR}/ankerl/include")

# using seqan
#target_include_directories (hixf PUBLIC 
#    ${CMAKE_CURRENT_SOURCE_DIR} 
#    ${CMAKE_CURRENT_BINARY_DIR}/../seqan/lemon/include 
#    ${CMAKE_CURRENT_BINARY_DIR}/../seqan/chopper/include
#    )
target_link_libraries (hashing unordered_dense::unordered_dense seqan3::seqan3)

#install(TARGETS IBF DESTINATION lib COMPONENT libraries)

#install(FILES ${IBF_HDRS} DESTINATION include COMPONENT headers)
#install(FILES ${WEIGHTS} DESTINATION data COMPONENT addfiles)
