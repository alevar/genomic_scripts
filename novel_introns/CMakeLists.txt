cmake_minimum_required(VERSION 3.10)
project(novel_introns)

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

include_directories(src)

add_executable(novel_introns
        src/arg_parse.cpp
        src/arg_parse.h
        src/codons.cpp
        src/codons.h
        src/GBase.cpp
        src/GBase.h
        src/gdna.cpp
        src/gdna.h
        src/GFaSeqGet.cpp
        src/GFaSeqGet.h
        src/GHash.hh
        src/GList.hh
        src/GVec.hh
        src/gff.h
        src/gff.cpp
        src/tokenize.cpp
        src/TrackingStats.cpp
        src/TrackingStats.h
        src/TrackingTree.cpp
        src/TrackingTree.h
        main.cpp)

target_link_libraries(novel_introns -lbz2 -llzma -pthread -lcrypto -lcurl -lz)

add_executable(novel_introns1
        src/arg_parse.cpp
        src/arg_parse.h
        src/codons.cpp
        src/codons.h
        src/GBase.cpp
        src/GBase.h
        src/gdna.cpp
        src/gdna.h
        src/GFaSeqGet.cpp
        src/GFaSeqGet.h
        src/GHash.hh
        src/GList.hh
        src/GVec.hh
        src/gff.h
        src/gff.cpp
        src/tokenize.cpp
        src/TrackingStats.cpp
        src/TrackingStats.h
        src/TrackingTree.cpp
        src/TrackingTree.h
        novel_introns1.cpp)

target_link_libraries(novel_introns1 -lbz2 -llzma -pthread -lcrypto -lcurl -lz)