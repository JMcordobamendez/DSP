cmake_minimum_required(VERSION 3.16)
project(pyfda_cpp VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(PYFDA_BUILD_GUI "Build the Qt GUI" ON)

if(MSVC)
    add_compile_options(/W3 /utf-8 /permissive-)
    add_compile_definitions(_USE_MATH_DEFINES)
else()
    add_compile_options(-Wall -Wextra)
endif()

# Numerical core without any Qt dependency
add_library(pyfda_core STATIC
    src/core/conversions.cpp
    src/core/data_io.cpp
    src/core/filter_design.cpp
    src/core/filter_io.cpp
    src/core/filtering.cpp
    src/core/fir_design.cpp
    src/core/iir_design.cpp
    src/core/poly.cpp
    src/core/remez.cpp
    src/core/special.cpp
    src/core/stimulus.cpp
)
target_include_directories(pyfda_core PUBLIC src/core)

# CLI used for the comparison with scipy (tests/verify_scipy.py)
add_executable(pyfda_cli tools/pyfda_cli.cpp)
target_link_libraries(pyfda_cli PRIVATE pyfda_core)

if(PYFDA_BUILD_GUI)
    set(CMAKE_AUTOMOC ON)
    find_package(Qt6 REQUIRED COMPONENTS Widgets)
    add_executable(pyfda_cpp WIN32
        src/gui/main.cpp
        src/gui/main_window.cpp
        src/gui/spec_panel.cpp
        src/gui/plot_widget.cpp
        src/gui/response_tabs.cpp
        src/gui/data_filt_tab.cpp
        src/gui/coeffs_tab.cpp
        src/gui/tran_tab.cpp
        src/gui/coeffs_tab.hpp
        src/gui/data_filt_tab.hpp
        src/gui/design_view.hpp
        src/gui/logger.hpp
        src/gui/main_window.hpp
        src/gui/plot_widget.hpp
        src/gui/response_tabs.hpp
        src/gui/spec_panel.hpp
        src/gui/tran_tab.hpp
    )
    target_link_libraries(pyfda_cpp PRIVATE pyfda_core Qt6::Widgets)
    install(TARGETS pyfda_cpp RUNTIME DESTINATION .)
endif()
