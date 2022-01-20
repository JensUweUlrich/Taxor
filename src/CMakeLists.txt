
cmake_minimum_required (VERSION 3.16)

project (PoreClassify CXX C)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)




if (WIN32)
    if (MSVC)
        add_compile_options(-D_WIN32_WINNT=0x600)
        add_compile_options("/MT")
        add_compile_options("/std:c++20")
        add_compile_options("/Zc:__cplusplus")
        add_compile_options("/Zc:__cpp_concepts")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++20 /MT /Zc:__cplusplus /Zc:__cpp_concepts")
    else()
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++17")
        add_compile_options("--std=c++17")
    endif()
#    set(CMAKE_ASM_NASM_OBJECT_FORMAT win64)
    set(CMAKE_INSTALL_PREFIX "C:/Apps/${PROJECT_NAME}")
    
else()
#    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
    set(CMAKE_INSTALL_PREFIX "/usr/local/${PROJECT_NAME}")
    add_compile_options("-Ofast")
#    add_link_options("stdc++fs")
#    add_compile_options("-std=c++17")
endif()

# add all subdirectories

add_subdirectory(seqan)
#add_subdirectory(boost)
add_subdirectory(main)

