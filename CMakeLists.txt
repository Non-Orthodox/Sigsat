find_package(Eigen3 3.3 REQUIRED NO_MODULE)
if(TARGET Eigen3::Eigen)
  message(STATUS "Eigen3 v${EIGEN3_VERSION_STRING} found in ${EIGEN3_INCLUDE_DIR}")
endif()

set(CORE src/gps_common.cpp
          src/binary_ops.cpp
          src/gps_ephemeris.cpp
          src/gps_lnav_data.cpp
          src/gps_correlator_sim.cpp
          src/gps_signal_gen.cpp
  )
add_library(Sigsat ${CORE})
target_include_directories(Sigsat
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
  )
target_link_libraries(Sigsat PUBLIC Eigen3::Eigen)

add_subdirectory(unit_tests)

