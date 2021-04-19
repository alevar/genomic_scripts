cmake_minimum_required(VERSION 3.10)
project(orfanage)

set (CMAKE_STATIC_LINKER_FLAGS "-Wl,--as-needed.-lcurl")

#option(TIEBRUSH_BUILD_LIBBIGWIG "Turn on/off building step of libBigWig (used for bioconda)." ON)
#if (NOT TIEBRUSH_BUILD_LIBBIGWIG)
#    add_definitions(-DTIEBRUSH_BUILD_LIBBIGWIG=0)
#    set(LIBBIGWIG_MAKE_CMD "")
#else ()
#    add_definitions(-DTIEBRUSH_BUILD_LIBBIGWIG=1)
#    set(LIBBIGWIG_MAKE_CMD "make")
#endif ()

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

include_directories(src)
include_directories("${CMAKE_SOURCE_DIR}/include/")

#ExternalProject_Add(libBigWig
#        GIT_REPOSITORY "https://github.com/alevar/libBigWig"
#        UPDATE_COMMAND ""
#        DOWNLOAD_DIR "${CMAKE_SOURCE_DIR}/include/"
#        SOURCE_DIR "${CMAKE_SOURCE_DIR}/include/libBigWig"
#        BUILD_IN_SOURCE 1
#        CONFIGURE_COMMAND ""
#        BUILD_COMMAND "${LIBBIGWIG_MAKE_CMD}"
#        INSTALL_COMMAND ""
#        )
#
#include_directories(${CMAKE_SOURCE_DIR}/include/libBigWig)
#link_directories(${CMAKE_SOURCE_DIR}/include/libBigWig)

add_executable(orfanage
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

add_dependencies(orfanage libBigWig)
#target_link_libraries(orfanage ${CMAKE_SOURCE_DIR}/include/libBigWig/libBigWig.a -lbz2 -llzma -pthread -lcrypto -lcurl -lz)

target_link_libraries(orfanage -lbz2 -llzma -pthread -lcrypto -lcurl -lz)