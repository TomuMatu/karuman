cmake_minimum_required(VERSION 3.13)

include(pico_sdk_import.cmake)

project(karman C CXX ASM)
set(CMAKE_C_STNDARD 11)
set(CMAKE_CXX_STANDARD 17)

pico_sdk_init()

include_directories(/home/tomumatu/pico/eigen)

add_executable(karman
        karman.cpp
)

# Pull in our pico_stdlib which pulls in commonly used features
target_link_libraries(karman pico_stdlib hardware_uart)

# enable usb output, disable uart output
pico_enable_stdio_usb(karman 1)
pico_enable_stdio_uart(karman 0)

# create map/bin/hex file etc.
pico_add_extra_outputs(karman)
