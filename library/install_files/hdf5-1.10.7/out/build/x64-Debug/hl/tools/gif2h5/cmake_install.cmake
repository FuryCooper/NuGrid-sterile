# Install script for directory: C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/hl/tools/gif2h5

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xhltoolsapplicationsx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/bin/gif2h5.exe")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xhltoolsapplicationsx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/bin/gif2h5-shared.exe")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xhltoolsapplicationsx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/bin/h52gif.exe")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xhltoolsapplicationsx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/bin/h52gif-shared.exe")
endif()

