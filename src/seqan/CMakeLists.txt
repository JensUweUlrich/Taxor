message(STATUS "${CMAKE_GENERATOR}")

if (WIN32)
  set(sdsl_cache_args  "-DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite
                        -DMSVC_COMPILER")
else()
  set(sdsl_cache_args  "-DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite")
endif()


configure_file(CMakeLists.txt.in ${CMAKE_CURRENT_BINARY_DIR}/seqan2-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/seqan2-download)
if(result)
  message(FATAL_ERROR "CMake step for seqan failed: ${result}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build . --config Release
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/seqan2-download )

add_library(sdsl STATIC IMPORTED GLOBAL)
if (WIN32)
  set_target_properties(sdsl PROPERTIES 
    IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/seqan2-download/sdsl-lite-prefix/src/sdsl-lite-build/lib/Release/sdsl.lib")
else()
  set_target_properties(sdsl PROPERTIES 
    IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/seqan2-download/sdsl-lite-prefix/src/sdsl-lite-build/lib/libsdsl.a")
endif()
target_include_directories( sdsl SYSTEM INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite/include )

if (WIN32)
set(ZLIB_ROOT "${CMAKE_CURRENT_BINARY_DIR}/zlib")
add_library(zlib STATIC IMPORTED GLOBAL)
if (WIN32)
    set_target_properties(zlib PROPERTIES IMPORTED_LOCATION "${ZLIB_ROOT}/lib/zlibstatic.lib")
    set(ZLIB_LIBRARIES zlib)
else()
    set_target_properties(zlib PROPERTIES IMPORTED_LOCATION "${ZLIB_ROOT}/lib/libz.a")
    set(ZLIB_LIBRARIES zlib)
endif()

set(ZLIB_INCLUDE_DIRS "${ZLIB_ROOT}/include")
endif()

find_package(ZLIB REQUIRED)

set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_BINARY_DIR}/seqan2/seqan2-src/util/cmake)
set(SEQAN_INCLUDE_PATH ${CMAKE_CURRENT_BINARY_DIR}/seqan2/seqan2-src/include)
set(SeqAn_DIR ${CMAKE_CURRENT_BINARY_DIR}/seqan2/seqan2-src/util/cmake)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

find_package( SeqAn CONFIG REQUIRED PATHS ${SeqAn_DIR} )

add_library( seqan INTERFACE )
target_include_directories( seqan SYSTEM INTERFACE ${SEQAN_INCLUDE_DIRS} )
target_link_libraries( seqan INTERFACE ${SEQAN_LIBRARIES} )

# Add include directories.

#message(FATAL_ERROR ${SEQAN_INCLUDE_DIRS})
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

