#ifndef SATELLITE_CONSTELLATIONS_INCLUDE_GPS_EPHEMERIS
#define SATELLITE_CONSTELLATIONS_INCLUDE_GPS_EPHEMERIS

#include <random>

#include <Eigen/Dense>

#include "common_types.hpp"

namespace Gps
{

uint32_t ParamToBinary(const double param, const double scale_factor);

double ParamFromBinary(const uint32_t data, const double scale_factor, const uint8_t num_bits,
  const bool twos_comp);


struct ClockData
{
  double T_GD {0.0};
  double t_oc {0.0};
  double a_f0 {0.0};
  double a_f1 {0.0};
  double a_f2 {0.0};
  uint16_t IODC {0};

  double Offset(const double gps_time) const;
  double OffsetRate(const double gps_time) const;
  double OffsetRateRate() const;

  void Randomize();
};


struct Ephemeris
{
  double M_0 {0.0};
  double del_n {0.0};
  double e {0.0};
  double sqrtA {0.0};
  double Omega_0 {0.0};
  double i_0 {0.0};
  double omega {0.0};
  double Omega_dot {0.0};
  double IDOT {0.0};
  
  double C_uc {0.0};
  double C_us {0.0};
  double C_rc {0.0};
  double C_rs {0.0};
  double C_ic {0.0};
  double C_is {0.0};

  double t_oe {0.0};
  uint8_t IODE {0};
  
  inline double EfromAnomaly(const double M_k, const unsigned int iterations) const;
  inline double EfromTime(const double gps_time, const unsigned int iterations) const;

  void P(const double gps_time, Eigen::Vector3d& pos) const;
  void PV(const double gps_time, Eigen::Vector3d& pos, Eigen::Vector3d& vel) const;
  void PVA(const double gps_time, Eigen::Vector3d& pos, Eigen::Vector3d& vel,
    Eigen::Vector3d& acc) const;
  void PVA(const double gps_time, Eigen::Vector3d& pos, Eigen::Vector3d& vel,
    Eigen::Vector3d& acc, bool calc_vel, bool calc_accel) const;

  double RelTime(const double gps_time) const;
  double RelTimeRate(const double gps_time) const;
  double RelTimeRateRate(const double gps_time) const;

  void Randomize();
  
  constexpr static double J2 = 0.0010826262;
  constexpr static double RELETIVISTIC_F = -4.442807633e-10;
  constexpr static double WGS84_MU = 3.986005e14;
  constexpr static double WGS84_EARTH_RATE = 7.2921151467e-5;
  constexpr static double WGS84_EQUAT_RADIUS = 6378137.0;
};

constexpr ClockData ClockDataScaleFactors =
{
  std::pow(2.0,-31),
  16,
  std::pow(2.0,-31),
  std::pow(2.0,-43),
  std::pow(2.0,-55),
  1
};

constexpr ClockData ClockDataLowerLimits = 
{
  -std::pow(2.0,7-31),
  0,
  -std::pow(2.0,21-31),
  -std::pow(2.0,15-43),
  -std::pow(2.0,7-55),
  0
};

constexpr ClockData ClockDataUpperLimits = 
{
  127 * std::pow(2.0,-31),
  604784,
  (std::pow(2.0,21) - 1.0) * std::pow(2.0,-31),
  (std::pow(2.0,15) - 1.0) * std::pow(2.0,-43),
  127 * std::pow(2.0,-55),
  1023
};

constexpr Ephemeris EphemerisScaleFactors =
{
  std::pow(2.0,-31),
  std::pow(2.0,-43),
  std::pow(2.0,-33),
  std::pow(2.0,-19),
  std::pow(2.0,-31),
  std::pow(2.0,-31),
  std::pow(2.0,-31),
  std::pow(2.0,-43),
  std::pow(2.0,-43),
  std::pow(2.0,-29),
  std::pow(2.0,-29),
  std::pow(2.0,-5),
  std::pow(2.0,-5),
  std::pow(2.0,-29),
  std::pow(2.0,-29),
  16,
  1
};

constexpr Ephemeris EphemerisLowerLimits = 
{
  -std::pow(2.0,31-31),
  -std::pow(2.0,15-43),
  0.0,
  2530.0,
  -std::pow(2.0,31-31),
  -std::pow(2.0,31-31),
  -std::pow(2.0,31-31),
  -6.33e-7,
  -std::pow(2.0,13-43), // IDOT
  -std::pow(2.0,15-29), // C_uc
  -std::pow(2.0,15-29), // C_us
  -std::pow(2.0,15-5), // C_rc
  -std::pow(2.0,15-5), // C_rs
  -std::pow(2.0,15-29), // C_ic
  -std::pow(2.0,15-29), // C_is
  0,
  0
};

constexpr Ephemeris EphemerisUpperLimits = 
{
  (std::pow(2.0,31) - 1.0) * std::pow(2.0,-31),
  (std::pow(2.0,15) - 1.0) * std::pow(2.0,-43),
  0.03,
  8192.0,
  (std::pow(2.0,31) - 1.0) * std::pow(2.0,-31),
  (std::pow(2.0,31) - 1.0) * std::pow(2.0,-31),
  (std::pow(2.0,31) - 1.0) * std::pow(2.0,-31),
  0.0,
  (std::pow(2.0,13) - 1.0) * std::pow(2.0,-43), // IDOT
  (std::pow(2.0,15) - 1.0) * std::pow(2.0,-29), // C_uc
  (std::pow(2.0,15) - 1.0) * std::pow(2.0,-29), // C_us
  (std::pow(2.0,15) - 1.0) * std::pow(2.0,-5), // C_rc
  (std::pow(2.0,15) - 1.0) * std::pow(2.0,-5), // C_rs
  (std::pow(2.0,15) - 1.0) * std::pow(2.0,-29), // C_ic
  (std::pow(2.0,15) - 1.0) * std::pow(2.0,-29), // C_is
  604784,
  255
};


} // namespace gps
#endif
