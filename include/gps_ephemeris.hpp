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

  double RelTime(const double gps_time) const;
  double RelTimeRate(const double gps_time) const;
  double RelTimeRateRate(const double gps_time) const;

  void Randomize();

  void Print() const;
  
  constexpr static double J2 = 0.0010826262;
  constexpr static double RELETIVISTIC_F = -4.442807633e-10;
  constexpr static double WGS84_MU = 3.986005e14;
  constexpr static double WGS84_EARTH_RATE = 7.2921151467e-5;
  constexpr static double WGS84_EQUAT_RADIUS = 6378137.0;

  // void PVA(const double gps_time, Eigen::Vector3d& pos, Eigen::Vector3d& vel,
  //   Eigen::Vector3d& acc, bool calc_vel, bool calc_accel) const;
  
  template<bool CalcVel, bool CalcAccel>
  void CalcPVA(const double gps_time, Eigen::Vector3d& pos, Eigen::Vector3d& vel,
    Eigen::Vector3d& accel) const
  {
    double A = std::pow(sqrtA,2.0);
    
    double n_0 = std::sqrt(WGS84_MU / std::pow(A,3.0));
    double t_k = gps_time - t_oe;
    if (t_k > 302400.0) {
      t_k -= 604800.0;
    } else if (t_k < -302400) {
      t_k += 604800.0;
    }
    // std::cout << "t_k: " << t_k << '\n';
    
    double n = n_0 + del_n;
    // std::cout << "n: " << n << '\n';
    double M_k = M_0 + (n * t_k);
    M_k = std::fmod(M_k + TwoPi<double>, TwoPi<double>);
    // std::cout << "M: " << M_k << '\n';

    double E_k = EfromAnomaly(M_k,5);
    E_k = std::fmod(E_k + TwoPi<double>, TwoPi<double>);
    // std::cout << "E_k: " << E_k << '\n';

    // double v_k = std::sqrt((1.0+e)/(1.0-e)) * std::tan(E_k / 2.0);
    // v_k = 2.0 * std::atan(v_k);
    double cv_k = (std::cos(E_k) - e) / (1.0 - (e * std::cos(E_k)));
    double sv_k = (std::sqrt(1.0 - (e*e)) * std::sin(E_k)) / (1.0 - (e * std::cos(E_k)));
    double v_k = std::atan2(sv_k,cv_k);
    // std::cout << "v_k: " << v_k << '\n';

    double Phi_k = v_k + omega;
    circular_fmod(Phi_k,TwoPi<double>);
    // std::cout << "Phi_k: " << Phi_k << '\n';

    double Phi2 = 2.0 * Phi_k;
    double du_k = (C_us * std::sin(Phi2)) + (C_uc * std::cos(Phi2));
    double dr_k = (C_rs * std::sin(Phi2)) + (C_rc * std::cos(Phi2));
    double di_k = (C_is * std::sin(Phi2)) + (C_ic * std::cos(Phi2));
    // std::cout << "du_k: " << du_k << '\n';
    // std::cout << "dr_k: " << dr_k << '\n';
    // std::cout << "di_k: " << di_k << '\n';
    
    double u_k = Phi_k + du_k;
    double r_k = A * (1.0 - (e * std::cos(E_k))) + dr_k;
    double i_k = i_0 + di_k + (IDOT * t_k);
    // std::cout << "u_k: " << u_k << '\n';
    // std::cout << "r_k: " << r_k << '\n';
    // std::cout << "i_k: " << i_k << '\n';

    double x_orb = r_k * std::cos(u_k);
    double y_orb = r_k * std::sin(u_k);

    double Omega_k = Omega_0 + ((Omega_dot - WGS84_EARTH_RATE) * t_k) - (WGS84_EARTH_RATE * t_oe);
    circular_fmod(Omega_k, TwoPi<double>);
    // std::cout << "Omega_k: " << Omega_k << '\n';

    pos(0) = (x_orb * std::cos(Omega_k)) - (y_orb * std::cos(i_k) * std::sin(Omega_k));
    pos(1) = (x_orb * std::sin(Omega_k)) + (y_orb * std::cos(i_k) * std::cos(Omega_k));
    pos(2) = y_orb * std::sin(i_k);

    // Velocity Terms
    if constexpr (CalcVel) {
      double Ed_k = n / (1.0 - (e * std::cos(E_k)));
      double vd_k = Ed_k * std::sqrt(1.0 - (e*e)) / (1.0 - (e * std::cos(E_k)));
      double id_k = IDOT + ( 2.0 * vd_k * ((C_is * std::cos(Phi2)) - (C_ic * std::sin(Phi2))) );
      double ud_k = vd_k + ( 2.0 * vd_k * ((C_us * std::cos(Phi2)) - (C_uc * std::sin(Phi2))) );
      double rd_k = (e * A * Ed_k * std::sin(E_k))
                    + ( 2.0 * vd_k * ((C_rs * std::cos(Phi2)) - (C_rc * std::sin(Phi2))) );
      double Omega_dot_k = Omega_dot - WGS84_EARTH_RATE;
      
      double xd_orb = (rd_k * std::cos(u_k)) - (r_k * ud_k * std::sin(u_k));
      double yd_orb = (rd_k * std::sin(u_k)) + (r_k * ud_k * std::cos(u_k));
      
      vel(0) = (-x_orb * Omega_dot_k * std::sin(Omega_k))
              + (xd_orb * std::cos(Omega_k))
              - (yd_orb * std::sin(Omega_k) * std::cos(i_k))
              - ( y_orb * ((Omega_dot_k * std::cos(Omega_k) * std::cos(i_k)) 
                        - (id_k * std::sin(Omega_k * std::sin(i_k)))) ); 
    
      vel(1) = (x_orb * Omega_dot_k * std::cos(Omega_k))
              + (xd_orb * std::sin(Omega_k))
              + (yd_orb * std::cos(Omega_k) * std::cos(i_k))
              - ( y_orb * ((Omega_dot_k * std::sin(Omega_k) * std::cos(i_k)) 
                        + (id_k * std::cos(Omega_k) * std::sin(i_k))) ); 
    
      vel(2) = (yd_orb * std::sin(i_k)) + (y_orb * id_k * std::cos(i_k));
      
    }

    // Acceleration Terms
    if constexpr (CalcAccel) {  
      double r2 = r_k * r_k;
      double r3 = r2 * r_k;
      double F = -1.5 * J2 * (WGS84_MU / r2) * std::pow(WGS84_EQUAT_RADIUS / r_k, 2.0);
      double F_term = F * ( 1.0 - (5.0 * std::pow(pos(2) / r_k, 2.0)) );  
      double omega_e2 = std::pow(WGS84_EARTH_RATE, 2.0);

      accel(0) = (-WGS84_MU * pos(0) / r3) + (F_term * pos(0) / r_k)
                + (2.0 * vel(1) * WGS84_EARTH_RATE) + (pos(0) * omega_e2);
    
      accel(1) = (-WGS84_MU * pos(1) / r3) + (F_term * pos(1) / r_k)
                - (2.0 * vel(0) * WGS84_EARTH_RATE) + (pos(1) * omega_e2);
    
      accel(2) = (-WGS84_MU * pos(2) / r3)
                + ( F * (3.0 - (5.0 * std::pow(pos(2) / r_k, 2.0))) * pos(2) / r_k );  
    }
  }
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
