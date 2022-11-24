add_custom_target(
    kepler_air
    COMMAND xcrun -sdk macosx metal -c ${CMAKE_CURRENT_SOURCE_DIR}/kepler.metal -o ${CMAKE_CURRENT_BINARY_DIR}/kepler.air
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/kepler.metal
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/kepler.air)
add_custom_target(
    kepler_metallib
    COMMAND xcrun -sdk macosx metallib ${CMAKE_CURRENT_BINARY_DIR}/kepler.air -o ${CMAKE_CURRENT_BINARY_DIR}/default.metallib
    DEPENDS kepler_air
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/default.metallib)
