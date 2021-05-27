# Search for header files in specific paths in home directory

if(NOT Blitz_INCLUDE_DIRS)
    set(_blitz_INCLUDE_SEARCH_DIRS "")

    if(BLITZ_INCLUDEDIR)
        list(APPEND _boost_INCLUDE_SEARCH_DIRS ${BLITZ_INCLUDEDIR})
    endif()

    if(BLITZ_ROOT)
        list(APPEND _blitz_INCLUDE_SEARCH_DIRS ${BLITZ_ROOT}/include ${BLITZ_ROOT})
        if(BLITZ_INCLUDEDIR)
            list(APPEND _blitz_INCLUDE_SEARCH_DIRS ${BLITZ_ROOT}/${BLITZ_INCLUDEDIR})
        endif()
    endif()

    list(APPEND _blitz_INCLUDE_SEARCH_DIRS $ENV{HOME} PATH_SUFFIXES /local/include /include /env/include)

    if( Blitz_NO_SYSTEM_PATHS)
        list(APPEND _blitz_INCLUDE_SEARCH_DIRS NO_CMAKE_SYSTEM_PATH NO_SYSTEM_ENVIRONMENT_PATH)
    endif()

    find_path(Blitz_INCLUDE_DIRS
        NAMES         blitz/blitz.h
        HINTS         ${_blitz_INCLUDE_SEARCH_DIRS}
    )
endif()

if(NOT Blitz_LIBRARIES)
    set(_blitz_LIBRARY_SEARCH_DIRS "")

    if(BLITZ_LIBRARYDIR)
        list(APPEND _boost_LIBRARY_SEARCH_DIRS ${BLITZ_LIBRARYDIR})
    endif()

    if(BLITZ_ROOT)
        list(APPEND _blitz_LIBRARY_SEARCH_DIRS ${BLITZ_ROOT}/lib ${BLITZ_ROOT}/lib64 ${BLITZ_ROOT})
        if(BLITZ_LIBRARYDIR)
            list(APPEND _blitz_LIBRARY_SEARCH_DIRS ${BLITZ_ROOT}/${BLITZ_LIBRARYDIR})
        endif()
    endif()

    list(APPEND _blitz_LIBRARY_SEARCH_DIRS $ENV{HOME} /lib /lib64 /local/lib /local/lib64 /env/lib /env/lib64) 

    if( Blitz_NO_SYSTEM_PATHS)
        list(APPEND _blitz_LIBRARY_SEARCH_DIRS NO_CMAKE_SYSTEM_PATH NO_SYSTEM_ENVIRONMENT_PATH)
    endif()

    find_library(Blitz_LIBRARIES
        NAMES blitz
        HINTS ${_blitz_LIBRARY_SEARCH_DIRS}
    )
endif()

if(Blitz_INCLUDE_DIRS)
    set(Blitz_INCLUDE_FOUND TRUE)
    set(Blitz_FOUND TRUE)
endif()

if(Blitz_LIBRARIES)
    set(Blitz_LIBRARIES_FOUND TRUE)
endif()

if(Blitz_INCLUDE_FOUND)
    if(NOT Blitz_FIND_QUIETLY)
        message(STATUS "Found Blitz headers: ${Blitz_INCLUDE_DIRS}")
    endif()
else()
    if(Blitz_FIND_REQUIRED)
        message(FATAL_ERROR "Could not find Blitz (which is required)")
    endif()
endif()

if((CMAKE_BUILD_TYPE MATCHES Debug) OR (CMAKE_BUILD_TYPE MATCHES PIGSDebug))
    if(Blitz_LIBRARIES_FOUND)
        if(NOT Blitz_FIND_QUIETLY)
            message(STATUS "Found Blitz libraries: ${Blitz_LIBRARIES}")
        endif()
    else()
        message(FATAL_ERROR "Could not find Blitz libraries (which are required for debug builds)")
    endif()
endif()
