cmake_minimum_required(VERSION 3.10)
project(ender)

set (CMAKE_STATIC_LINKER_FLAGS "-Wl,--as-needed.-lcurl")

include(ExternalProject)

# thanks for the macro: David Grayson
macro(use_cxx11)
    if (CMAKE_VERSION VERSION_LESS "3.1")
        if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
        endif ()
    else ()
        set (CMAKE_CXX_STANDARD 11)
    endif ()
endmacro(use_cxx11)

#use_cxx11()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
set (CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fno-omit-frame-pointer -fsanitize=address")


include_directories(src)
include_directories("${CMAKE_SOURCE_DIR}/include/")

add_executable(ender
        src/arg_parse.cpp
        src/arg_parse.h
        ${CMAKE_SOURCE_DIR}/include/gclib/gff.h
        ${CMAKE_SOURCE_DIR}/include/gclib/gff.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/codons.h
        ${CMAKE_SOURCE_DIR}/include/gclib/codons.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/gdna.h
        ${CMAKE_SOURCE_DIR}/include/gclib/gdna.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GVec.hh
        ${CMAKE_SOURCE_DIR}/include/gclib/GList.hh
        ${CMAKE_SOURCE_DIR}/include/gclib/GHash.hh
        ${CMAKE_SOURCE_DIR}/include/gclib/GFaSeqGet.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GFaSeqGet.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GFastaIndex.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GFastaIndex.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GFastaFile.h
        ${CMAKE_SOURCE_DIR}
        src/gff_utils.h
        src/gff_utils.cpp
        main.cpp)

target_link_libraries(ender -lbz2 -llzma -pthread -lcrypto -lcurl -lz)
