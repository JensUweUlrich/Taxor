set(Main_SRC "main.cpp")
file(GLOB MAIN_HDRS "*.hpp")
add_executable (${PROJECT_NAME} ${MAIN_HDRS} ${Main_SRC})
#find_package (Threads)


# using seqan3
if (WIN32)
	add_definitions(-DMSVC_COMPILER)
	target_link_libraries (${PROJECT_NAME} sdsl seqan)
else()
	target_link_libraries (${PROJECT_NAME} sdsl seqan)
endif()
