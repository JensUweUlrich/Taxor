
find_package(zlib REQUIRED)
if (WIN32)
    set(BOOST_CONFIGURE_CMD "bootstrap.bat")
    set(BOOST_BUILD_FLAGS "toolset=msvc 
                    link=static 
                    variant=release 
                    threading=multi 
                    runtime-link=static
                    debug-symbols=off
                    optimization=speed 
                    cxxflags=\"/MT\"
                    address-model=64
                    -s ZLIB_INCLUDE=${CMAKE_CURRENT_BINARY_DIR}/zlib/include
                    -s ZLIB_LIBPATH==${CMAKE_CURRENT_BINARY_DIR}/zlib/lib
                    -s ZLIB_SOURCE=${CMAKE_CURRENT_BINARY_DIR}/zlib
                    -s ZLIB_NAME=zlib
                    --prefix=${CMAKE_CURRENT_BINARY_DIR}/boost/build/install")
else()
    set(BOOST_CONFIGURE_CMD "bootstrap.sh")
    set(BOOST_BUILD_FLAGS "toolset=gcc 
                    link=static 
                    variant=release 
                    threading=multi 
                    runtime-link=static
                    debug-symbols=off
                    optimization=speed 
                    address-model=64
                    --prefix=${CMAKE_CURRENT_BINARY_DIR}/boost/build/install")
endif()

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

configure_file(CMakeLists.txt.in boost-build/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/boost-build)
if(result)
  message(FATAL_ERROR "CMake step for boost failed: ${result}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build . --config Release
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/boost-build)
if(result)
  message(FATAL_ERROR "Build step for boost failed: ${result}")
endif()

find_package(zlib REQUIRED)

set(BOOST_MIN_VERSION 1.54)
set(BOOST_ROOT "${CMAKE_CURRENT_BINARY_DIR}/boost")
set(Boost_USE_STATIC_LIBS ON)
set(BOOST_LIBRARY_DIRS "${BOOST_ROOT}/build/install/lib")

add_definitions( -DBOOST_ALL_NO_LIB )

if (WIN32)
	file(GLOB boost_libs "${BOOST_LIBRARY_DIRS}/*-s-x64-1_76.lib")
	set(Boost_INCLUDE_DIRS "${BOOST_ROOT}/build/install/include/boost-1_76")
else()
	file(GLOB boost_libs "${BOOST_LIBRARY_DIRS}/*.a")
	set(Boost_INCLUDE_DIRS "${BOOST_ROOT}/build/install/include")
endif()

include_directories(${Boost_INCLUDE_DIRS})

foreach(file ${boost_libs})
    set(Boost_LIBRARIES ${Boost_LIBRARIES} ${file})
endforeach()


set(BOOST_COMPONENTS
    date_time
    filesystem
    iostreams
    log
    log_setup
    program_options
    thread
)

foreach (component ${BOOST_COMPONENTS})
    add_library(boost_${component} STATIC IMPORTED GLOBAL)
    if (WIN32)
        set_target_properties(boost_${component} PROPERTIES 
	     IMPORTED_LOCATION "${BOOST_LIBRARY_DIRS}/libboost_${component}-vc142-mt-s-x64-1_76.lib")
    else()
	set_target_properties(boost_${component} PROPERTIES 
             IMPORTED_LOCATION "${BOOST_LIBRARY_DIRS}/libboost_${component}.a")
    endif()
endforeach()

find_package(zlib REQUIRED)
message(STATUS ${Boost_INCLUDE_DIRS})