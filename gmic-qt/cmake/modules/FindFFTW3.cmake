# - Try to find the Fftw3 Libraries
#
# Once done this will define
#
#  FFTW3_FOUND        - system has fftw3
#  FFTW3_INCLUDE_DIRS - the fftw3 include directories
#  FFTW3_LIBRARIES    - the libraries needed to use fftw3
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

if (NOT WIN32)

    include(LibFindMacros)
    libfind_pkg_check_modules(FFTW3_PKGCONF fftw3>=3.2)

    find_path(FFTW3_INCLUDE_DIR
        NAMES fftw3.h
        HINTS ${FFTW3_PKGCONF_INCLUDE_DIRS} ${FFTW3_PKGCONF_INCLUDEDIR}
        PATH_SUFFIXES fftw3
    )

    find_library(FFTW3_LIBRARY_CORE
        NAMES fftw3
        HINTS ${FFTW3_PKGCONF_LIBRARY_DIRS} ${FFTW3_PKGCONF_LIBDIR}
    )

    set(FFTW3_CORE_PROCESS_INCLUDES FFTW3_INCLUDE_DIR)
    set(FFTW3_CORE_PROCESS_LIBS FFTW3_LIBRARY_CORE)
    libfind_process(FFTW3_CORE)

    if(FFTW3_CORE_FOUND)
        message(STATUS "FFTW core Found Libraries: " ${FFTW3_CORE_LIBRARIES})
    endif()

    find_library(FFTW3_LIBRARY_THREADS
        NAMES fftw3_threads
        HINTS ${FFTW3_PKGCONF_LIBRARY_DIRS} ${FFTW3_PKGCONF_LIBDIR}
    )

    set(FFTW3_THREADS_PROCESS_INCLUDES FFTW3_INCLUDE_DIR)
    set(FFTW3_THREADS_PROCESS_LIBS FFTW3_LIBRARY_THREADS)
    libfind_process(FFTW3_THREADS)

    if(FFTW3_THREADS_FOUND)
        message(STATUS "FFTW threads Found Libraries: " ${FFTW3_THREADS_LIBRARIES})
    endif()

    if(FFTW3_CORE_FOUND AND FFTW3_THREADS_FOUND)
        set(FFTW3_FOUND true)
        set(FFTW3_LIBRARIES ${FFTW3_CORE_LIBRARIES} ${FFTW3_THREADS_LIBRARIES})
    endif()

else()

    # TODO: Maybe use fftw3/FFTW3Config.cmake?

    find_path(FFTW3_INCLUDE_DIR
        NAMES fftw3.h
    )

    find_library(
        FFTW3_LIBRARY
        NAMES libfftw3 libfftw3-3 libfftw3f-3 libfftw3l-3 fftw3
        DOC "Libraries to link against for FFT Support")

    if (FFTW3_LIBRARY)
        set(FFTW3_LIBRARY_DIR ${FFTW3_LIBRARY})
    endif()

    set (FFTW3_LIBRARIES ${FFTW3_LIBRARY})

    if(FFTW3_INCLUDE_DIR AND FFTW3_LIBRARY_DIR)
       set (FFTW3_FOUND true)
       message(STATUS "Correctly found FFTW3")
    else()
       message(STATUS "Could not find FFTW3")
    endif()

endif()
