# Install script for directory: C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/install/x64-Debug")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/hdf5.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5api_adpt.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5public.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5Apublic.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5ACpublic.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5Cpublic.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5Dpublic.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5Epubgen.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5Epublic.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5Fpublic.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5FDcore.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5FDdirect.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5FDfamily.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5FDhdfs.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5FDlog.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5FDmirror.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5FDmpi.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5FDmpio.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5FDmulti.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5FDpublic.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5FDros3.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5FDs3comms.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5FDsec2.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5FDsplitter.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5FDstdio.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5FDwindows.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5Gpublic.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5Ipublic.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5Lpublic.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5MMpublic.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5Opublic.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5Ppublic.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5PLextern.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5PLpublic.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5Rpublic.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5Spublic.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5Tpublic.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5Zpublic.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5Epubgen.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5version.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/src/H5overflow.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/H5pubconf.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibrariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg]|[Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE FILE OPTIONAL FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/bin/.pdb")
  endif("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg]|[Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibrariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg]|[Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE FILE OPTIONAL FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/bin/.pdb")
  endif("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg]|[Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibrariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/bin/libhdf5_D.lib")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibrariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY OPTIONAL FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/bin/hdf5_D.lib")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibrariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE SHARED_LIBRARY FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/bin/hdf5_D.dll")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibrariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/CMakeFiles/hdf5-1.10.7.pc")
endif()

