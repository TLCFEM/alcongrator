cmake_minimum_required(VERSION 3.22)

project(alcongrator)

set(QT_PATH "$ENV{HOME}/Qt/6.9.0/gcc_64/lib/cmake" CACHE STRING "Path to Qt")
set(CMAKE_PREFIX_PATH "${QT_PATH}")
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOUIC_SEARCH_PATHS form)
set(CMAKE_CXX_STANDARD 17)

add_compile_definitions(ARMA_DONT_USE_ATLAS)
add_compile_definitions(ARMA_USE_SUPERLU)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")

if (UNIX)
    link_directories(lib/linux)
elseif (WIN32)
    if (MSVC)
        link_directories(lib/win-msvc)
    else ()
        link_directories(lib/win-gcc)
    endif ()
endif ()

include_directories(include)
include_directories(include/QCustomPlot)

find_package(Qt6 REQUIRED COMPONENTS
        Core
        Gui
        Widgets
        PrintSupport
        Multimedia
)

set(SOURCES
        include/QCustomPlot/qcustomplot.cpp
        src/alcongrator.cpp
        src/CustomFilter.cpp
        src/CubicSpline.cpp
        src/IntegrationScheme.cpp
        src/MainWindow.cpp
        src/resampling.cpp
)

set(UIS
        form/MainWindow.ui
)

set(RESOURCES
        res/res.qrc
)

include_directories(src)

add_subdirectory(superlu-src)

add_executable(${PROJECT_NAME} ${UIS} ${SOURCES} ${RESOURCES})

add_dependencies(${PROJECT_NAME} superlu)

target_link_libraries(${PROJECT_NAME}
        Qt6::Core
        Qt6::Gui
        Qt6::Widgets
        Qt6::PrintSupport
        Qt6::Multimedia
        superlu
        openblas
)

if (WIN32)
    target_link_libraries(${PROJECT_NAME} opengl32)
elseif (UNIX)
    target_link_libraries(${PROJECT_NAME} GL gfortran quadmath gomp)
endif ()
