# Install script for directory: C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xcppheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5AbstractDs.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5Alltypes.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5ArrayType.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5AtomType.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5Attribute.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5Classes.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5CommonFG.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5CompType.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5Cpp.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5CppDoc.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5DataSet.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5DataSpace.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5DataType.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5DaccProp.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5DcreatProp.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5DxferProp.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5EnumType.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5Exception.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5FaccProp.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5FcreatProp.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5File.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5FloatType.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5Group.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5IdComponent.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5Include.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5IntType.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5LaccProp.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5LcreatProp.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5Library.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5Location.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5Object.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5OcreatProp.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5PredType.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5PropList.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5StrType.h"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/c++/src/H5VarLenType.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xcpplibrariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg]|[Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE FILE OPTIONAL FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/bin/.pdb")
  endif("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg]|[Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xcpplibrariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg]|[Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE FILE OPTIONAL FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/bin/.pdb")
  endif("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg]|[Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xcpplibrariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/bin/libhdf5_cpp_D.lib")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xcpplibrariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY OPTIONAL FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/bin/hdf5_cpp_D.lib")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xcpplibrariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE SHARED_LIBRARY FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/bin/hdf5_cpp_D.dll")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xcpplibrariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/CMakeFiles/hdf5_cpp-1.10.7.pc")
endif()

