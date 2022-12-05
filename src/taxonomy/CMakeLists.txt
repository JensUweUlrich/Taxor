file(GLOB TAXON_SRC "*.cpp" "*.c")
file(GLOB TAXON_HDRS "*.hpp" "*.h")

add_library (taxonomy STATIC ${TAXON_SRC} ${TAXON_HDRS})

# include spd fast logging
#set(spdlog_DIR ${CMAKE_CURRENT_BINARY_DIR}/../spdlog/spdlog/lib/cmake/spdlog)
#find_package(spdlog REQUIRED)

# include SafeQueue and StopClock
#include_directories("${CMAKE_CURRENT_BINARY_DIR}/../seqan/lemon/include" "${CMAKE_CURRENT_BINARY_DIR}/../seqan/chopper/include")

# using seqan
#target_include_directories (hixf PUBLIC 
#    ${CMAKE_CURRENT_SOURCE_DIR} 
#    ${CMAKE_CURRENT_BINARY_DIR}/../seqan/lemon/include 
#    ${CMAKE_CURRENT_BINARY_DIR}/../seqan/chopper/include
#    )
target_link_libraries (taxonomy seqan3::seqan3)

#install(TARGETS IBF DESTINATION lib COMPONENT libraries)

#install(FILES ${IBF_HDRS} DESTINATION include COMPONENT headers)
#install(FILES ${WEIGHTS} DESTINATION data COMPONENT addfiles)
