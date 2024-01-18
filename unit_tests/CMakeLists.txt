add_executable(gps_ca_tests gps_ca_tests.cpp)
target_link_libraries(gps_ca_tests PUBLIC Sigsat Eigen3::Eigen python_plotting)

