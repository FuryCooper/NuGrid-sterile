# Install script for directory: C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xconfiginstallx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/hdf5/hdf5-targets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/hdf5/hdf5-targets.cmake"
         "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/CMakeFiles/Export/cmake/hdf5/hdf5-targets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/hdf5/hdf5-targets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/hdf5/hdf5-targets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/hdf5" TYPE FILE FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/CMakeFiles/Export/cmake/hdf5/hdf5-targets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/hdf5" TYPE FILE FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/CMakeFiles/Export/cmake/hdf5/hdf5-targets-debug.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xconfiginstallx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/hdf5" TYPE FILE FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/CMakeFiles/hdf5-config.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xconfiginstallx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/hdf5" TYPE FILE FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/CMakeFiles/hdf5-config-version.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibrariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE FILE FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/libhdf5.settings")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xhdfdocumentsx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/." TYPE FILE FILES "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/COPYING")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xhdfdocumentsx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/." TYPE FILE FILES
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/release_docs/USING_HDF5_CMake.txt"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/release_docs/COPYING"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/release_docs/RELEASE.txt"
    "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/release_docs/USING_HDF5_VS.txt"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/src/cmake_install.cmake")
  include("C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/test/cmake_install.cmake")
  include("C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/utils/cmake_install.cmake")
  include("C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/tools/cmake_install.cmake")
  include("C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/examples/cmake_install.cmake")
  include("C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/hl/cmake_install.cmake")
  include("C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/c++/cmake_install.cmake")
  include("C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/hl/c++/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "C:/Users/13973/Desktop/NuGrid/NuGrid/library/install_files/hdf5-1.10.7/out/build/x64-Debug/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
