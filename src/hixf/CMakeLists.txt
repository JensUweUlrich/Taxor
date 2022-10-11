file(GLOB HIXF_SRC "*.cpp")
file(GLOB HIXF_HDRS "*.hpp")

add_library (hixf STATIC ${HIXF_SRC} ${HIXF_HDRS})

# include spd fast logging
#set(spdlog_DIR ${CMAKE_CURRENT_BINARY_DIR}/../spdlog/spdlog/lib/cmake/spdlog)
#find_package(spdlog REQUIRED)

# include SafeQueue and StopClock
include_directories("${CMAKE_CURRENT_BINARY_DIR}/../seqan/lemon/include" 
                    "${CMAKE_CURRENT_BINARY_DIR}/../seqan/chopper/include"
                    )

# using seqan
target_include_directories (hixf PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR} 
    ${CMAKE_CURRENT_BINARY_DIR}/../seqan/lemon/include 
    ${CMAKE_CURRENT_BINARY_DIR}/../seqan/chopper/include
    ${CMAKE_CURRENT_SOURCE_DIR}/../hashing
    )
target_link_libraries (hixf seqan3::seqan3 hashing)

#install(TARGETS IBF DESTINATION lib COMPONENT libraries)

#install(FILES ${IBF_HDRS} DESTINATION include COMPONENT headers)
#install(FILES ${WEIGHTS} DESTINATION data COMPONENT addfiles)
