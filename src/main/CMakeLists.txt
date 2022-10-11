file(GLOB Main_SRC "*.cpp" "*.c")
file(GLOB MAIN_HDRS "*.hpp" "*.h")
add_executable (${PROJECT_NAME} ${MAIN_HDRS} ${Main_SRC})
#find_package (Threads)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../hashing")

target_include_directories(${PROJECT_NAME} INTERFACE 
	${CMAKE_CURRENT_BINARY_DIR}/../hashing
	${CMAKE_CURRENT_BINARY_DIR}/../hixf
	)

# using seqan3
if (WIN32)
	#add_definitions(-DMSVC_COMPILER)
	target_link_libraries (${PROJECT_NAME} seqan3::seqan3)
else()
	target_link_libraries (${PROJECT_NAME} seqan3::seqan3 hashing hixf)
endif()
