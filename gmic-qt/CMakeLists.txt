project(gmic-qt)

message(STATUS "Using CMake version: ${CMAKE_VERSION}")

cmake_minimum_required(VERSION 3.10.0 FATAL_ERROR)
LIST (APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
include(FeatureSummary)
include(FindPkgConfig)

set(CMAKE_CXX_STANDARD 11)
add_definitions(-Dcimg_use_cpp11=1)
set(MIN_QT_VERSION 5.2.0)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC OFF)
set(CMAKE_AUTORCC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(EXTRA_LIBRARIES)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

message("Build type is " ${CMAKE_BUILD_TYPE})

set (GMIC_QT_HOST "gimp" CACHE STRING "Define for which host qmic-qt will be built: gimp, gimp3 (experimental), krita-plugin, none, paintdotnet or 8bf.")
if (${GMIC_QT_HOST} STREQUAL "none")
  message("Building standalone version.")
else()
  message("Building for target host application: " ${GMIC_QT_HOST})
endif()

option(ENABLE_SYSTEM_GMIC "Find GMIC shared library installed on the system" ON)

if (ENABLE_SYSTEM_GMIC)
  option(ENABLE_DYNAMIC_LINKING "Dynamically link the binaries to the GMIC shared library" ON)
else()
  option(ENABLE_DYNAMIC_LINKING "Dynamically link the binaries to the GMIC shared library" OFF)

  set (GMIC_LIB_PATH "${GMIC_PATH}" CACHE STRING "Define the path to the GMIC shared library")
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../src/gmic.cpp")
  set (GMIC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../src" CACHE STRING "Define the path to the gmic headers")
  else()
  set (GMIC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../gmic/src" CACHE STRING "Define the path to the gmic headers")
  endif()

  message("G'MIC path: " ${GMIC_PATH})
endif()

option(ENABLE_CURL "Add support for curl" ON)

option(ENABLE_ASAN "Enable -fsanitize=address (if debug build)" ON)
option(ENABLE_FFTW3 "Enable FFTW3 library support" ON)

include(CheckIPOSupported)

check_ipo_supported(RESULT HAVE_LTO LANGUAGES CXX)

if (HAVE_LTO)
option(ENABLE_LTO "Enable Link Time Optimizer" OFF)
endif()

if (MSVC)
  option(ENABLE_CFG "Enable Control Flow Guard (MSVC)" ON)
  add_definitions(-D__PRETTY_FUNCTION__=__FUNCSIG__)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(ENABLE_LTO OFF FORCE)
endif()


if (NOT ENABLE_SYSTEM_GMIC)
  #
  # Look for G'MIC repository
  #
  get_filename_component(GMIC_ABSOLUTE_PATH ${GMIC_PATH} ABSOLUTE BASEDIR ${CMAKE_SOURCE_DIR})
  if (EXISTS ${GMIC_ABSOLUTE_PATH}/gmic.cpp)
    message("Found G'MIC repository")
  else()
    get_filename_component(TARGET_CLONE_DIR ${GMIC_ABSOLUTE_PATH}/.. ABSOLUTE)
    message("")
    message("Cannot find G'MIC repository in " ${GMIC_ABSOLUTE_PATH} )
    message("")
    message("You should try:")
    message("")
    message("   git clone https://github.com/dtschump/gmic.git " ${TARGET_CLONE_DIR}/gmic )
    message("")
    message(FATAL_ERROR "\nG'MIC repository not found")
  endif()

  #
  # Look for CImg.h and gmic_stdlib_community.h
  #
  set(GMIC_FILES CImg.h gmic_stdlib_community.h)
  foreach(F ${GMIC_FILES})
  if(EXISTS ${GMIC_ABSOLUTE_PATH}/${F})
      message("Found " ${GMIC_PATH}/${F})
  else()
      message(${F} " not found")
      execute_process(COMMAND make -C ${GMIC_ABSOLUTE_PATH} ${F})
      if(EXISTS ${GMIC_ABSOLUTE_PATH}/${F})
      message("Found " ${GMIC_PATH}/${F})
      else()
      message(FATAL_ERROR "\nCannot obtain " ${GMIC_PATH}/${F})
      endif()
  endif()
  endforeach()

  #
  # Ensure that gmic and CImg are the same version
  #
  file(STRINGS ${GMIC_ABSOLUTE_PATH}/CImg.h CIMG_VERSION REGEX "cimg_version ")
  string(REGEX REPLACE ".*cimg_version " "" CIMG_VERSION ${CIMG_VERSION})
  message("CImg version is [" ${CIMG_VERSION} "]")

  file(STRINGS ${GMIC_ABSOLUTE_PATH}/gmic.h GMIC_VERSION REGEX "gmic_version ")
  string(REGEX REPLACE ".*gmic_version " "" GMIC_VERSION ${GMIC_VERSION})
  message("G'MIC version is [" ${GMIC_VERSION} "]")

  if (NOT(${GMIC_VERSION} EQUAL ${CIMG_VERSION}))
  message(FATAL_ERROR "\nVersion numbers of files 'gmic.h' (" ${GMIC_VERSION} ") and 'CImg.h' (" ${CIMG_VERSION} ") mismatch")
  endif()
endif()


option(PRERELEASE "Set to ON makes this a prelease build")
if (${PRERELEASE})
    string(TIMESTAMP PRERELEASE_DATE %y%m%d)
    message("Prelease date is " ${PRERELEASE_DATE})
    add_definitions(-Dgmic_prerelease="${PRERELEASE_DATE}")
endif()

option(DRMINGW "Set to ON enables the drmingw debugger.")
if (${DRMINGW})
    add_definitions(-DDRMINGW)
endif()


# Required packages

#
# Gmic
#
if (ENABLE_SYSTEM_GMIC)
  find_package(Gmic REQUIRED CONFIG)
endif (ENABLE_SYSTEM_GMIC)

#
# Threads
#
find_package(Threads REQUIRED)

#
# Qt
#
find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Core)
find_package(Qt${QT_VERSION_MAJOR} ${MIN_QT_VERSION}
        REQUIRED COMPONENTS
        Core
        Gui
        Widgets
        Network
        LinguistTools
)
if(QT_VERSION VERSION_GREATER_EQUAL "6.3")
  if(QT_NO_CREATE_VERSIONLESS_FUNCTIONS)
    qt6_standard_project_setup()
  else()
    qt_standard_project_setup()
  endif()
endif()

#
# PNG
#
find_package(PNG REQUIRED)
add_definitions(${PNG_DEFINITIONS})
add_definitions(-Dcimg_use_png)
include_directories(SYSTEM ${PNG_INCLUDE_DIR})
if (APPLE)
    # this is not added correctly on OSX -- see http://forum.kde.org/viewtopic.php?f=139&t=101867&p=221242#p221242
    include_directories(SYSTEM ${PNG_INCLUDE_DIR})
endif()

#
# ZLIB
#
find_package(ZLIB REQUIRED)
add_definitions(-Dcimg_use_zlib)
include_directories(SYSTEM ${ZLIB_INCLUDE_DIRS} )

#
# FFTW3
#
if (ENABLE_FFTW3)
  find_package(FFTW3 REQUIRED)
  add_definitions(-Dcimg_use_fftw3 )
  include_directories(${FFTW3_INCLUDE_DIR})

  # Detect
  include(CheckCXXSourceCompiles)
  set(CMAKE_REQUIRED_INCLUDES ${FFTW3_INCLUDE_DIR})
  set(CMAKE_REQUIRED_LIBRARIES ${FFTW3_LIBRARIES})
  check_cxx_source_compiles("
      #include <fftw3.h>

      int main() {
          fftw_init_threads();
      }
  " HAVE_FFTW3_THREADS)

  if(HAVE_FFTW3_THREADS)
    message(STATUS "FFTW threads Found")
    list(APPEND EXTRA_LIBRARIES ${FFTW3_THREADS_LIBRARIES})
  else()
    add_definitions(-Dcimg_use_fftw3_singlethread)
  endif()
endif()

#
# CURL
#
if(ENABLE_CURL)
    find_package(CURL)
    if (CURL_FOUND)
        add_definitions(-Dcimg_use_curl)
        include_directories(SYSTEM ${CURL_INCLUDE_DIRS} )
    endif()
endif()

#
# Test for OpenMP
#
find_package(OpenMP 2.0)
set_package_properties(OpenMP PROPERTIES
    DESCRIPTION "A low-level parallel execution library"
    URL "http://openmp.org/wp/"
    TYPE OPTIONAL
    PURPOSE "Optionally used by gmic-qt")

if (OpenMP_FOUND)
    message(STATUS "G'Mic: using OpenMP ${OpenMP_CXX_VERSION}")
    link_libraries(OpenMP::OpenMP_CXX)
    add_definitions(-Dcimg_use_openmp)
endif()

#
# LTO option
#

if (ENABLE_LTO)
    message(STATUS "Link Time Optimizer enabled")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

#
# Enable CFG
#
if (MSVC)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:16777216")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /STACK:16777216")
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} /STACK:16777216")
    if (ENABLE_CFG)
        add_compile_options(/guard:CF)
        add_link_options(/GUARD:CF)
    endif()
elseif(WIN32)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--stack,16777216")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--stack,16777216")
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Wl,--stack,16777216")
endif()

#
# add all defines
#

if (QT_NO_CREATE_VERSIONLESS_FUNCTIONS)
  set(gmic_qt_LIBRARIES
          Qt${QT_VERSION_MAJOR}::Core
          Qt${QT_VERSION_MAJOR}::Widgets
          Qt${QT_VERSION_MAJOR}::Gui
          Qt${QT_VERSION_MAJOR}::Network
  )
else()
  set(gmic_qt_LIBRARIES
          Qt::Core
          Qt::Widgets
          Qt::Gui
          Qt::Network
  )
endif()

set(gmic_qt_LIBRARIES
        ${gmic_qt_LIBRARIES}
        ${PNG_LIBRARIES}
        ${FFTW3_LIBRARIES}
        ${ZLIB_LIBRARIES}
        ${EXTRA_LIBRARIES}
)

if(ENABLE_CURL)
    if (CURL_FOUND)
        set(gmic_qt_LIBRARIES
            ${gmic_qt_LIBRARIES}
            ${CURL_LIBRARIES}
        )
    endif()
endif()

add_definitions(-Dgmic_core)
add_definitions(-Dgmic_community)
add_definitions(-Dcimg_use_abort)
add_definitions(-Dgmic_is_parallel)
add_definitions(-Dgmic_gui)
add_definitions(-Dcimg_use_abort)
add_definitions(-Dcimg_appname=\"gmic\")

if (UNIX)
    if(ANDROID)
        add_definitions(-Dcimg_display=0)
        add_definitions(-D_IS_UNIX_)
    elseif(NOT APPLE)
        add_definitions(-Dcimg_display=1)
        add_definitions(-D_IS_LINUX_)
        add_definitions(-Dcimg_use_vt100)
        add_definitions(-D_IS_UNIX_)
        find_package(X11)
        set(gmic_qt_LIBRARIES
            ${gmic_qt_LIBRARIES}
            ${X11_LIBRARIES} # XXX: Search for X11: Wayland is coming!
        )
    endif()
endif()

if (APPLE)
    add_definitions(-Dcimg_display=0)
    add_definitions(-D_IS_MACOS_)
    set(CMAKE_MACOSX_RPATH 1)
    set(BUILD_WITH_INSTALL_RPATH 1)
    add_definitions(-Wno-macro-redefined -Wno-deprecated-register)
    if (QT_VERSION_MAJOR VERSION_EQUAL "5")
      add_definitions(-mmacosx-version-min=10.9)
    endif()
endif()

if (WIN32)
    add_definitions(-Dcimg_display=2)
    add_definitions(-DPSAPI_VERSION=1)
    add_definitions(-D_IS_WINDOWS_)
    if (MSVC)
      add_definitions(-D_CRT_SECURE_NO_WARNINGS)
      add_compile_options(/wd4267)
    endif()
    set(gmic_qt_LIBRARIES
        ${gmic_qt_LIBRARIES}
        Threads::Threads psapi gdi32
    )
endif()

SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
SET(CMAKE_INSTALL_RPATH "$ORIGIN/")

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-D_GMIC_QT_DEBUG_)
    if(ENABLE_ASAN)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
    endif(ENABLE_ASAN)
elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
    add_definitions(-DQT_NO_DEBUG_OUTPUT)
    if (MSVC)
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /fp:fast /Oi")
    else()
      string(REPLACE "-O2" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
      string(REPLACE "-O3" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")
    endif()
    if (NOT MSVC)
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -s")
    endif()
    if (WIN32 AND NOT MSVC)
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mwindows")
    endif()
elseif (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    add_definitions(-DQT_NO_DEBUG_OUTPUT)
    if(MSVC)
      string(REPLACE "Ob1" "Ob2" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} /fp:fast /Oi")
    else()
      string(REPLACE "-O2" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
      string(REPLACE "-O3" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
      if (NOT ENABLE_SYSTEM_GMIC)
          set_source_files_properties(${GMIC_PATH}/gmic.cpp PROPERTIES COMPILE_FLAGS "-Ofast")
      endif (NOT ENABLE_SYSTEM_GMIC)
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O2")
    endif()
else()
    message(FATAL_ERROR "Build type not recognized (${CMAKE_BUILD_TYPE})")
endif()

if (NOT ENABLE_SYSTEM_GMIC)
    include_directories(${GMIC_PATH})
endif (NOT ENABLE_SYSTEM_GMIC)
include_directories(${CMAKE_SOURCE_DIR}/src)

set (gmic_qt_SRCS
  src/ClickableLabel.h
  src/Common.h
  src/CroppedActiveLayerProxy.h
  src/CroppedImageListProxy.h
  src/DialogSettings.h
  src/FilterParameters/AbstractParameter.h
  src/FilterParameters/BoolParameter.h
  src/FilterParameters/ButtonParameter.h
  src/FilterParameters/ChoiceParameter.h
  src/FilterParameters/ColorParameter.h
  src/FilterParameters/ConstParameter.h
  src/FilterParameters/CustomDoubleSpinBox.h
  src/FilterParameters/CustomSpinBox.h
  src/FilterParameters/FileParameter.h
  src/FilterParameters/FilterParametersWidget.h
  src/FilterParameters/FloatParameter.h
  src/FilterParameters/FolderParameter.h
  src/FilterParameters/IntParameter.h
  src/FilterParameters/LinkParameter.h
  src/FilterParameters/MultilineTextParameterWidget.h
  src/FilterParameters/NoteParameter.h
  src/FilterParameters/PointParameter.h
  src/FilterParameters/SeparatorParameter.h
  src/FilterParameters/TextParameter.h
  src/FilterSelector/FavesModel.h
  src/FilterSelector/FavesModelReader.h
  src/FilterSelector/FavesModelWriter.h
  src/FilterSelector/FiltersModelBinaryReader.h
  src/FilterSelector/FiltersModelBinaryWriter.h
  src/FilterSelector/FiltersModel.h
  src/FilterSelector/FiltersModelReader.h
  src/FilterSelector/FiltersPresenter.h
  src/FilterSelector/FiltersView/FiltersView.h
  src/FilterSelector/FiltersView/FilterTreeAbstractItem.h
  src/FilterSelector/FiltersView/FilterTreeFolder.h
  src/FilterSelector/FiltersView/FilterTreeItemDelegate.h
  src/FilterSelector/FiltersView/FilterTreeItem.h
  src/FilterSelector/FiltersView/TreeView.h
  src/FilterSelector/FiltersVisibilityMap.h
  src/FilterSelector/FilterTagMap.h
  src/FilterSyncRunner.h
  src/FilterTextTranslator.h
  src/FilterThread.h
  src/Globals.h
  src/GmicProcessor.h
  src/GmicQt.h
  src/GmicStdlib.h
  src/HeadlessProcessor.h
  src/Host/GmicQtHost.h
  src/HtmlTranslator.h
  src/IconLoader.h
  src/ImageTools.h
  src/InputOutputState.h
  src/KeypointList.h
  src/LanguageSettings.h
  src/LayersExtentProxy.h
  src/Logger.h
  src/MainWindow.h
  src/Misc.h
  src/OverrideCursor.h
  src/ParametersCache.h
  src/PersistentMemory.h
  src/Settings.h
  src/SourcesWidget.h
  src/Tags.h
  src/TimeLogger.h
  src/Updater.h
  src/Utils.h
  src/Widgets/InOutPanel.h
  src/Widgets/LanguageSelectionWidget.h
  src/Widgets/PreviewWidget.h
  src/Widgets/ProgressInfoWidget.h
  src/Widgets/ProgressInfoWindow.h
  src/Widgets/SearchFieldWidget.h
  src/Widgets/VisibleTagSelector.h
  src/Widgets/ZoomLevelSelector.h
  src/ZoomConstraint.h
)

if (NOT ENABLE_SYSTEM_GMIC)
    set(gmic_qt_SRCS
        ${gmic_qt_SRCS}
        ${GMIC_PATH}/gmic.h
        ${GMIC_PATH}/CImg.h
        ${GMIC_PATH}/gmic_stdlib_community.h
    )
endif()

set(gmic_qt_SRCS
  ${gmic_qt_SRCS}
  src/ClickableLabel.cpp
  src/Common.cpp
  src/CroppedActiveLayerProxy.cpp
  src/CroppedImageListProxy.cpp
  src/DialogSettings.cpp
  src/FilterParameters/AbstractParameter.cpp
  src/FilterParameters/BoolParameter.cpp
  src/FilterParameters/ButtonParameter.cpp
  src/FilterParameters/ChoiceParameter.cpp
  src/FilterParameters/ColorParameter.cpp
  src/FilterParameters/ConstParameter.cpp
  src/FilterParameters/CustomDoubleSpinBox.cpp
  src/FilterParameters/CustomSpinBox.cpp
  src/FilterParameters/FileParameter.cpp
  src/FilterParameters/FilterParametersWidget.cpp
  src/FilterParameters/FloatParameter.cpp
  src/FilterParameters/FolderParameter.cpp
  src/FilterParameters/IntParameter.cpp
  src/FilterParameters/LinkParameter.cpp
  src/FilterParameters/MultilineTextParameterWidget.cpp
  src/FilterParameters/NoteParameter.cpp
  src/FilterParameters/PointParameter.cpp
  src/FilterParameters/SeparatorParameter.cpp
  src/FilterParameters/TextParameter.cpp
  src/FilterSelector/FavesModel.cpp
  src/FilterSelector/FavesModelReader.cpp
  src/FilterSelector/FavesModelWriter.cpp
  src/FilterSelector/FiltersModelBinaryReader.cpp
  src/FilterSelector/FiltersModelBinaryWriter.cpp
  src/FilterSelector/FiltersModel.cpp
  src/FilterSelector/FiltersModelReader.cpp
  src/FilterSelector/FiltersPresenter.cpp
  src/FilterSelector/FiltersView/FiltersView.cpp
  src/FilterSelector/FiltersView/FilterTreeAbstractItem.cpp
  src/FilterSelector/FiltersView/FilterTreeFolder.cpp
  src/FilterSelector/FiltersView/FilterTreeItem.cpp
  src/FilterSelector/FiltersView/FilterTreeItemDelegate.cpp
  src/FilterSelector/FiltersView/TreeView.cpp
  src/FilterSelector/FiltersVisibilityMap.cpp
  src/FilterSelector/FilterTagMap.cpp
  src/FilterSyncRunner.cpp
  src/FilterTextTranslator.cpp
  src/FilterThread.cpp
  src/Globals.cpp
  src/GmicProcessor.cpp
  src/GmicQt.cpp
  src/GmicStdlib.cpp
  src/HeadlessProcessor.cpp
  src/HtmlTranslator.cpp
  src/IconLoader.cpp
  src/ImageTools.cpp
  src/InputOutputState.cpp
  src/KeypointList.cpp
  src/LanguageSettings.cpp
  src/LayersExtentProxy.cpp
  src/Logger.cpp
  src/MainWindow.cpp
  src/Misc.cpp
  src/OverrideCursor.cpp
  src/ParametersCache.cpp
  src/PersistentMemory.cpp
  src/Settings.cpp
  src/SourcesWidget.cpp
  src/Tags.cpp
  src/TimeLogger.cpp
  src/Updater.cpp
  src/Utils.cpp
  src/Widgets/InOutPanel.cpp
  src/Widgets/LanguageSelectionWidget.cpp
  src/Widgets/PreviewWidget.cpp
  src/Widgets/ProgressInfoWidget.cpp
  src/Widgets/ProgressInfoWindow.cpp
  src/Widgets/SearchFieldWidget.cpp
  src/Widgets/VisibleTagSelector.cpp
  src/Widgets/ZoomLevelSelector.cpp
)

set (gmic_qt_FORMS
  ui/dialogsettings.ui
  ui/filtersview.ui
  ui/headlessprogressdialog.ui
  ui/inoutpanel.ui
  ui/languageselectionwidget.ui
  ui/mainwindow.ui
  ui/multilinetextparameterwidget.ui
  ui/progressinfowidget.ui
  ui/progressinfowindow.ui
  ui/SearchFieldWidget.ui
  ui/sourceswidget.ui
  ui/zoomlevelselector.ui
)

if(ENABLE_DYNAMIC_LINKING)
  set(CMAKE_SKIP_RPATH TRUE)
  # G'MIC-Qt needs visibility into the private symbols defined
  # by the gmic.cpp plugin. However, this is only possible
  # if the library is static OR if it's dynamic and built by
  # a compiler that supports .so-style exports.
  if (TARGET libgmicstatic OR MSVC OR NOT ENABLE_SYSTEM_GMIC)
    set(gmic_qt_LIBRARIES
      ${gmic_qt_LIBRARIES}
      libgmicstatic
    )
  elseif(TARGET libgmic)
    set(gmic_qt_LIBRARIES
      ${gmic_qt_LIBRARIES}
      libgmic
    )
  elseif(GMIC_LIB_PATH)
    set(gmic_qt_LIBRARIES
      ${gmic_qt_LIBRARIES}
      "gmic"
    )
  else()
    message(FATAL_ERROR "No G'MIC library is available for linking. Please build libgmic as a static library.")
  endif()
  if (NOT ENABLE_SYSTEM_GMIC)
    if (GMIC_LIB_PATH)
      link_directories(${GMIC_LIB_PATH})
      # Inject the G'MIC CImg plugin.
      include_directories(../src)
    else()
      # Mimic an external G'MIC library build for catching link ABI errors.
      add_library(libgmicstatic STATIC ../src/gmic.cpp)
      target_include_directories(libgmicstatic PUBLIC ../src)
      # We need internal access into the gmic-core API.
      target_compile_definitions(libgmicstatic PUBLIC gmic_core)
      set_target_properties(libgmicstatic
        PROPERTIES
          AUTOMOC OFF
      )
      target_link_libraries(libgmicstatic PUBLIC
        ${PNG_LIBRARIES}
        ${FFTW3_LIBRARIES}
        ${ZLIB_LIBRARIES}
        ${CURL_LIBRARIES}
        ${EXTRA_LIBRARIES})
    endif()
  else()
    # Inject the G'MIC CImg plugin.
    include_directories(../src)
  endif()
else(ENABLE_DYNAMIC_LINKING)
  set(gmic_qt_SRCS
    ${gmic_qt_SRCS}
    ${GMIC_PATH}/gmic.cpp
    )
endif(ENABLE_DYNAMIC_LINKING)

if (UNIX OR APPLE OR MSYS)
  message("Producing translation .qm files")
  execute_process(COMMAND make WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/translations OUTPUT_QUIET)
  message("Producing filter translation .qm files")
  execute_process(COMMAND make WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/translations/filters OUTPUT_QUIET)
else()
  set(gmic_translation_files
      translations/cs.ts
      translations/de.ts
      translations/es.ts
      translations/fr.ts
      translations/id.ts
      translations/it.ts
      translations/ja.ts
      translations/nl.ts
      translations/pl.ts
      translations/pt.ts
      translations/ru.ts
      translations/sv.ts
      translations/uk.ts
      translations/zh.ts
      translations/zh_tw.ts
  )

  set_source_files_properties(${gmic_translation_files} PROPERTIES OUTPUT_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/translations)

  set(gmic_filters_translation_files
      translations/filters/de.ts
      translations/filters/es.ts
      translations/filters/fr.ts
      translations/filters/it.ts
      translations/filters/ja.ts
      translations/filters/nl.ts
      translations/filters/pl.ts
      translations/filters/pt.ts
      translations/filters/ru.ts
      translations/filters/zh.ts
  )

  set_source_files_properties(${gmic_filters_translation_files} PROPERTIES OUTPUT_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/translations/filters)

  if(QT_VERSION VERSION_LESS_EQUAL "5.15.0")
    qt5_add_translation(gmic_qt_QM
      ${gmic_translation_files}
      ${gmic_filters_translation_files}
    )
  elseif(QT_NO_CREATE_VERSIONLESS_FUNCTIONS)
    qt6_add_translation(gmic_qt_QM
      ${gmic_translation_files}
      ${gmic_filters_translation_files}
    )
  else()
    qt_add_translation(gmic_qt_QM
      ${gmic_translation_files}
      ${gmic_filters_translation_files}
    )
  endif()
endif()

include(GNUInstallDirs)

install(FILES ${CMAKE_SOURCE_DIR}/../resources/gmic_cluts.gmz DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/gmic)
install(FILES ${CMAKE_SOURCE_DIR}/../resources/gmic_denoise_cnn.gmz DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/gmic)

set(gmic_qt_QRC
    gmic_qt.qrc
    translations.qrc
    wip_translations.qrc
)

if (${GMIC_QT_HOST} STREQUAL "gimp" OR ${GMIC_QT_HOST} STREQUAL "gimp3")
    if (${GMIC_QT_HOST} STREQUAL "gimp3")
        set(TARGET_GIMP_VERSION 3)
    else()
        set(TARGET_GIMP_VERSION 2)
    endif()

    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GIMP REQUIRED gimp-${TARGET_GIMP_VERSION}.0 IMPORTED_TARGET)
    # CMake does not support passing --define-variable through pkg_get_variable.
    execute_process(COMMAND ${PKG_CONFIG_EXECUTABLE} gimp-${TARGET_GIMP_VERSION}.0 --define-variable=prefix=${CMAKE_INSTALL_PREFIX} --variable gimplibdir OUTPUT_VARIABLE GIMP_PKGLIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)

    set (gmic_qt_SRCS ${gmic_qt_SRCS} src/Host/Gimp/host_gimp.cpp)
    if(QT_VERSION VERSION_LESS_EQUAL "5.15.0")
      qt5_wrap_ui(gmic_qt_SRCS ${gmic_qt_FORMS})
    elseif(QT_NO_CREATE_VERSIONLESS_FUNCTIONS)
      qt6_wrap_ui(gmic_qt_SRCS ${gmic_qt_FORMS})
    else()
      qt_wrap_ui(gmic_qt_SRCS ${gmic_qt_FORMS})
    endif()
    add_definitions(-DGMIC_HOST=gimp -DGIMP_DISABLE_DEPRECATED)
    add_executable(gmic_gimp_qt ${gmic_qt_SRCS} ${gmic_qt_QRC} ${gmic_qt_QM})
    target_link_libraries(
      gmic_gimp_qt
      PRIVATE
      PkgConfig::GIMP
      ${gmic_qt_LIBRARIES}
      )
    install(TARGETS gmic_gimp_qt RUNTIME DESTINATION "${GIMP_PKGLIBDIR}/plug-ins/gmic_gimp_qt")

elseif (${GMIC_QT_HOST} STREQUAL "krita-plugin")
    set(MIN_FRAMEWORKS_VERSION 5.44.0)

    find_package(ECM 5.22 REQUIRED NOMODULE)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${ECM_MODULE_PATH} ${ECM_KDE_MODULE_DIR})

    include(KDEInstallDirs)
    include(KDECMakeSettings)

    if (ANDROID)
        set (KRITA_PLUGIN_INSTALL_DIR ${LIB_INSTALL_DIR})
    else()
        set (KRITA_PLUGIN_INSTALL_DIR ${LIB_INSTALL_DIR}/kritaplugins)
    endif()

    find_package(KF5 ${MIN_FRAMEWORKS_VERSION} REQUIRED COMPONENTS
        CoreAddons
    )

    message(STATUS "Looking for Krita QMic libraries in: ${CMAKE_PREFIX_PATH}")

    find_library(KIS_IMAGE_INTERFACE_LIBRARY
      NAMES kritaqmicinterface
      REQUIRED)

    find_path(KIS_IMAGE_INTERFACE_DIR
      NAMES kis_qmic_plugin_interface.h
      REQUIRED)

    set_package_properties(kritaqmicinterface PROPERTIES
                           URL "http://www.krita.org"
                           DESCRIPTION "Krita GMic core library"
    )

    set (gmic_qt_SRCS ${gmic_qt_SRCS} src/Host/KritaPlugin/host.cpp src/Host/KritaPlugin/gmicqttoolplugin.cpp)
    set (gmic_qt_SRCS ${gmic_qt_SRCS} )
    if(QT_VERSION VERSION_LESS_EQUAL "5.15.0")
      qt5_wrap_ui(gmic_qt_SRCS ${gmic_qt_FORMS})
    elseif(QT_NO_CREATE_VERSIONLESS_FUNCTIONS)
      qt6_wrap_ui(gmic_qt_SRCS ${gmic_qt_FORMS})
    else()
      qt_wrap_ui(gmic_qt_SRCS ${gmic_qt_FORMS})
    endif()
    add_definitions(-DGMIC_HOST=krita-plugin)
    add_definitions(-D_GMIC_QT_DISABLE_THEMING_)
    add_definitions(-D_GMIC_QT_CONSENT_TO_UPDATE_FIRST_)
    add_definitions(-D_GMIC_QT_DISABLE_TRANSLATION_)
    add_definitions(-D_GMIC_USE_HOSTED_SETTINGS_)
    add_definitions(-D_GMIC_ASYNC_IMAGE_UPDATE_)
    add_definitions(-D_GMIC_HOSTED_PLUGIN_)
    add_library(krita_gmic_qt MODULE ${gmic_qt_SRCS} ${gmic_qt_QRC} ${gmic_qt_QM})
    target_include_directories(
      krita_gmic_qt
      PUBLIC
      ${KIS_IMAGE_INTERFACE_DIR}
    )
    target_link_libraries(
      krita_gmic_qt
      PRIVATE
      ${gmic_qt_LIBRARIES}
      ${KIS_IMAGE_INTERFACE_LIBRARY}
      KF5::CoreAddons
      )
    if (ANDROID)
      target_link_libraries(
        krita_gmic_qt
        PRIVATE
        log
      )
    endif()
    install(TARGETS krita_gmic_qt DESTINATION ${KRITA_PLUGIN_INSTALL_DIR}) # plugin

elseif (${GMIC_QT_HOST} STREQUAL "none")

  set (gmic_qt_SRCS ${gmic_qt_SRCS}
    src/Host/None/host_none.cpp
    src/Host/None/ImageDialog.h
    src/Host/None/ImageDialog.cpp
    src/Host/None/JpegQualityDialog.h
    src/Host/None/JpegQualityDialog.cpp
    )
  set(gmic_qt_FORMS ${gmic_qt_FORMS}
    src/Host/None/jpegqualitydialog.ui
    )
  if(QT_VERSION VERSION_LESS "5.15.0")
    qt5_wrap_ui(gmic_qt_SRCS ${gmic_qt_FORMS})
  elseif(QT_NO_CREATE_VERSIONLESS_FUNCTIONS)
    qt6_wrap_ui(gmic_qt_SRCS ${gmic_qt_FORMS})
  else()
    qt_wrap_ui(gmic_qt_SRCS ${gmic_qt_FORMS})
  endif()
  add_definitions(-DGMIC_HOST=standalone)
  add_executable(gmic_qt ${gmic_qt_SRCS} ${gmic_qt_QRC} ${gmic_qt_QM})
  target_link_libraries(gmic_qt PRIVATE ${gmic_qt_LIBRARIES})
  install(TARGETS gmic_qt RUNTIME DESTINATION bin)

elseif (${GMIC_QT_HOST} STREQUAL "paintdotnet")

  set (gmic_qt_SRCS ${gmic_qt_SRCS} src/Host/PaintDotNet/host_paintdotnet.cpp)
  if(QT_VERSION VERSION_LESS "5.15.0")
    qt5_wrap_ui(gmic_qt_SRCS ${gmic_qt_FORMS})
  elseif(QT_NO_CREATE_VERSIONLESS_FUNCTIONS)
    qt6_wrap_ui(gmic_qt_SRCS ${gmic_qt_FORMS})
  else()
    qt_wrap_ui(gmic_qt_SRCS ${gmic_qt_FORMS})
  endif()
  add_definitions(-DGMIC_HOST=paintdotnet)
  add_executable(gmic_paintdotnet_qt ${gmic_qt_SRCS} ${gmic_qt_QRC} ${gmic_qt_QM})
  target_link_libraries(
    gmic_paintdotnet_qt
    PRIVATE ${gmic_qt_LIBRARIES}
    )

elseif (${GMIC_QT_HOST} STREQUAL "8bf")

  # Look for a CMake package on MSVC or a PkgConfig file on MinGW etc.
  if (MSVC)
    find_package(lcms2 CONFIG REQUIRED)
    include_directories(${LCMS2_INCLUDE_DIR})
  else()
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(LCMS2 REQUIRED lcms2)
  endif()

  set (gmic_qt_SRCS ${gmic_qt_SRCS} src/Host/8bf/host_8bf.cpp)
  if(QT_VERSION VERSION_LESS_EQUAL "5.15.0")
    qt5_wrap_ui(gmic_qt_SRCS ${gmic_qt_FORMS})
  elseif(QT_NO_CREATE_VERSIONLESS_FUNCTIONS)
    qt6_wrap_ui(gmic_qt_SRCS ${gmic_qt_FORMS})
  else()
    qt_wrap_ui(gmic_qt_SRCS ${gmic_qt_FORMS})
  endif()
  add_definitions(-DGMIC_HOST=plugin8bf)
  add_executable(gmic_8bf_qt ${gmic_qt_SRCS} ${gmic_qt_QRC} ${gmic_qt_QM})
  target_link_libraries(
    gmic_8bf_qt
    PRIVATE
    ${gmic_qt_LIBRARIES}
	${LCMS2_LIBRARIES}
    )

else()
    message(FATAL_ERROR "GMIC_QT_HOST is not defined as gimp, gimp3, kritaplugin, none, paintdotnet or 8bf")
endif()

feature_summary(WHAT ALL FATAL_ON_MISSING_REQUIRED_PACKAGES)
