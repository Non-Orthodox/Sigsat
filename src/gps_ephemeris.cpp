#include <cmath>
#include <numbers>
#include <random>
#include <iostream>

#include <Eigen/Dense>

#include "gps_common.hpp"
#include "gps_ephemeris.hpp"
#include "binary_ops.hpp"

namespace Gps
{

uint32_t ParamToBinary(const double param, const double scale_factor)
{
  double quotient = param / scale_factor;
  return static_cast<uint32_t>(quotient < 0 ? quotient - 0.5 : quotient + 0.5);
}


double ParamFromBinary(uint32_t data, const double scale_factor, const uint8_t num_bits,
  const bool twos_comp)
{
  assert((num_bits >= 1) && (num_bits <= 32));
  assert(scale_factor > 0.0);

  static const uint32_t one = 1; // just avoiding the differing default size of integer literals

  if (num_bits != 32)
    data &= ((one << num_bits) - 1); // sets all MSBs past (num_bits-1) position to zero
  if (twos_comp && (bitVal<true>(data, num_bits-1))) {
    uint32_t signed_bit = (one << (num_bits-1));
    data &= (~signed_bit);
    data = (signed_bit - data); // un-signs the integer
    return -scale_factor * static_cast<double>(data);
  } else {
    return scale_factor * static_cast<double>(data);
  }
}


// --------------------------- ClockData --------------------------- 
double ClockData::Offset(const double gps_time) const
{
  double dt = gps_time - t_oc;
  return a_f0 + (a_f1 * dt) + (a_f2 * dt * dt);
}

double ClockData::OffsetRate(const double gps_time) const
{
  return a_f1 + (2.0 * (gps_time - t_oc) * a_f2);
}

double ClockData::OffsetRateRate() const
{
  return 2.0 * a_f2;
}

void ClockData::Randomize()
{
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  T_GD = distribution(internal::random_gen) * (ClockDataUpperLimits.T_GD - ClockDataLowerLimits.T_GD) + ClockDataLowerLimits.T_GD;
  t_oc = distribution(internal::random_gen) * (ClockDataUpperLimits.t_oc - ClockDataLowerLimits.t_oc) + ClockDataLowerLimits.t_oc;
  a_f0 = distribution(internal::random_gen) * (ClockDataUpperLimits.a_f0 - ClockDataLowerLimits.a_f0) + ClockDataLowerLimits.a_f0;
  a_f1 = distribution(internal::random_gen) * (ClockDataUpperLimits.a_f1 - ClockDataLowerLimits.a_f1) + ClockDataLowerLimits.a_f1;
  a_f2 = distribution(internal::random_gen) * (ClockDataUpperLimits.a_f2 - ClockDataLowerLimits.a_f2) + ClockDataLowerLimits.a_f2;
}


// --------------------------- Ephemeris --------------------------- 
inline double Ephemeris::EfromAnomaly(const double M_k, const unsigned int iterations) const
{
  double E_k = M_k;
  // for (unsigned int j = 0; j < iterations; j++) {
  //   E_k += (M_k - E_k + (e * std::sin(E_k))) / (1.0 - (e * std::cos(E_k)));
  // }

  double del_E;
  do {
    del_E = (M_k - E_k + (e * std::sin(E_k))) / (1.0 - (e * std::cos(E_k)));
    E_k += del_E;
  }
  while (del_E > 1.0e-15);
  return E_k;
}

inline double Ephemeris::EfromTime(const double gps_time, const unsigned int iterations) const
{
  double A = std::pow(sqrtA,2.0);
  
  double n_0 = std::sqrt(WGS84_MU / std::pow(A,3.0));
  double t_k = gps_time - t_oe;
  if (t_k > 302400.0) {
    t_k -= 604800.0;
  } else if (t_k < -302400) {
    t_k += 604800.0;
  }
  
  double n = n_0 + del_n;
  return EfromAnomaly(M_0 + (n * t_k),5);
}

void Ephemeris::P(const double gps_time, Eigen::Vector3d& pos) const
{
  Eigen::Vector3d filler;
  CalcPVA<false,false>(gps_time, pos, filler, filler);
}

void Ephemeris::PV(const double gps_time, Eigen::Vector3d& pos, Eigen::Vector3d& vel) const
{
  Eigen::Vector3d filler;
  CalcPVA<true,false>(gps_time, pos, vel, filler);
}

void Ephemeris::PVA(const double gps_time, Eigen::Vector3d& pos, Eigen::Vector3d& vel,
  Eigen::Vector3d& accel) const
{
  CalcPVA<true,true>(gps_time, pos, vel, accel);
}

// //TODO template<bool CalcVel, bool CalcAccel>
// void Ephemeris::PVA(const double gps_time, Eigen::Vector3d& pos, Eigen::Vector3d& vel,
//   Eigen::Vector3d& accel, bool calc_vel, bool calc_accel) const
// {
//   double A = std::pow(sqrtA,2.0);
  
//   double n_0 = std::sqrt(WGS84_MU / std::pow(A,3.0));
//   double t_k = gps_time - t_oe;
//   if (t_k > 302400.0) {
//     t_k -= 604800.0;
//   } else if (t_k < -302400) {
//     t_k += 604800.0;
//   }
//   // std::cout << "t_k: " << t_k << '\n';
  
//   double n = n_0 + del_n;
//   // std::cout << "n: " << n << '\n';
//   double M_k = M_0 + (n * t_k);
//   M_k = std::fmod(M_k + TwoPi<double>, TwoPi<double>);
//   // std::cout << "M: " << M_k << '\n';

//   double E_k = EfromAnomaly(M_k,5);
//   E_k = std::fmod(E_k + TwoPi<double>, TwoPi<double>);
//   // std::cout << "E_k: " << E_k << '\n';

//   // double v_k = std::sqrt((1.0+e)/(1.0-e)) * std::tan(E_k / 2.0);
//   // v_k = 2.0 * std::atan(v_k);
//   double cv_k = (std::cos(E_k) - e) / (1.0 - (e * std::cos(E_k)));
//   double sv_k = (std::sqrt(1.0 - (e*e)) * std::sin(E_k)) / (1.0 - (e * std::cos(E_k)));
//   double v_k = std::atan2(sv_k,cv_k);
//   // std::cout << "v_k: " << v_k << '\n';

//   double Phi_k = v_k + omega;
//   circular_fmod(Phi_k,TwoPi<double>);
//   // std::cout << "Phi_k: " << Phi_k << '\n';

//   double Phi2 = 2.0 * Phi_k;
//   double du_k = (C_us * std::sin(Phi2)) + (C_uc * std::cos(Phi2));
//   double dr_k = (C_rs * std::sin(Phi2)) + (C_rc * std::cos(Phi2));
//   double di_k = (C_is * std::sin(Phi2)) + (C_ic * std::cos(Phi2));
//   // std::cout << "du_k: " << du_k << '\n';
//   // std::cout << "dr_k: " << dr_k << '\n';
//   // std::cout << "di_k: " << di_k << '\n';
  
//   double u_k = Phi_k + du_k;
//   double r_k = A * (1.0 - (e * std::cos(E_k))) + dr_k;
//   double i_k = i_0 + di_k + (IDOT * t_k);
//   // std::cout << "u_k: " << u_k << '\n';
//   // std::cout << "r_k: " << r_k << '\n';
//   // std::cout << "i_k: " << i_k << '\n';

//   double x_orb = r_k * std::cos(u_k);
//   double y_orb = r_k * std::sin(u_k);

//   double Omega_k = Omega_0 + ((Omega_dot - WGS84_EARTH_RATE) * t_k) - (WGS84_EARTH_RATE * t_oe);
//   circular_fmod(Omega_k, TwoPi<double>);
//   // std::cout << "Omega_k: " << Omega_k << '\n';

//   pos(0) = (x_orb * std::cos(Omega_k)) - (y_orb * std::cos(i_k) * std::sin(Omega_k));
//   pos(1) = (x_orb * std::sin(Omega_k)) + (y_orb * std::cos(i_k) * std::cos(Omega_k));
//   pos(2) = y_orb * std::sin(i_k);

//   // Velocity Terms
//   if (calc_vel) {
//     double Ed_k = n / (1.0 - (e * std::cos(E_k)));
//     double vd_k = Ed_k * std::sqrt(1.0 - (e*e)) / (1.0 - (e * std::cos(E_k)));
//     double id_k = IDOT + ( 2.0 * vd_k * ((C_is * std::cos(Phi2)) - (C_ic * std::sin(Phi2))) );
//     double ud_k = vd_k + ( 2.0 * vd_k * ((C_us * std::cos(Phi2)) - (C_uc * std::sin(Phi2))) );
//     double rd_k = (e * A * Ed_k * std::sin(E_k))
//                   + ( 2.0 * vd_k * ((C_rs * std::cos(Phi2)) - (C_rc * std::sin(Phi2))) );
//     double Omega_dot_k = Omega_dot - WGS84_EARTH_RATE;
    
//     double xd_orb = (rd_k * std::cos(u_k)) - (r_k * ud_k * std::sin(u_k));
//     double yd_orb = (rd_k * std::sin(u_k)) + (r_k * ud_k * std::cos(u_k));
    
//     vel(0) = (-x_orb * Omega_dot_k * std::sin(Omega_k))
//             + (xd_orb * std::cos(Omega_k))
//             - (yd_orb * std::sin(Omega_k) * std::cos(i_k))
//             - ( y_orb * ((Omega_dot_k * std::cos(Omega_k) * std::cos(i_k)) 
//                       - (id_k * std::sin(Omega_k * std::sin(i_k)))) ); 
  
//     vel(1) = (x_orb * Omega_dot_k * std::cos(Omega_k))
//             + (xd_orb * std::sin(Omega_k))
//             + (yd_orb * std::cos(Omega_k) * std::cos(i_k))
//             - ( y_orb * ((Omega_dot_k * std::sin(Omega_k) * std::cos(i_k)) 
//                       + (id_k * std::cos(Omega_k) * std::sin(i_k))) ); 
  
//     vel(2) = (yd_orb * std::sin(i_k)) + (y_orb * id_k * std::cos(i_k));
    
//   }

//   // Acceleration Terms
//   if (calc_accel) {  
//     double r2 = r_k * r_k;
//     double r3 = r2 * r_k;
//     double F = -1.5 * J2 * (WGS84_MU / r2) * std::pow(WGS84_EQUAT_RADIUS / r_k, 2.0);
//     double F_term = F * ( 1.0 - (5.0 * std::pow(pos(2) / r_k, 2.0)) );  
//     double omega_e2 = std::pow(WGS84_EARTH_RATE, 2.0);

//     accel(0) = (-WGS84_MU * pos(0) / r3) + (F_term * pos(0) / r_k)
//               + (2.0 * vel(1) * WGS84_EARTH_RATE) + (pos(0) * omega_e2);
  
//     accel(1) = (-WGS84_MU * pos(1) / r3) + (F_term * pos(1) / r_k)
//               - (2.0 * vel(0) * WGS84_EARTH_RATE) + (pos(1) * omega_e2);
  
//     accel(2) = (-WGS84_MU * pos(2) / r3)
//               + ( F * (3.0 - (5.0 * std::pow(pos(2) / r_k, 2.0))) * pos(2) / r_k );  
//   }
// }


double Ephemeris::RelTime(const double gps_time) const
{
  return RELETIVISTIC_F * std::numbers::e * sqrtA * std::sin(EfromTime(gps_time,5));
}

double Ephemeris::RelTimeRate(const double gps_time) const
{
  double A = std::pow(sqrtA,2.0);
  double n_0 = std::sqrt(WGS84_MU / std::pow(A,3.0));
  double t_k = gps_time - t_oe;
  if (t_k > 302400.0) {
    t_k -= 604800.0;
  } else if (t_k < -302400) {
    t_k += 604800.0;
  }  
  double n = n_0 + del_n;
  double e_cos_E = std::numbers::e * std::cos( EfromAnomaly(M_0 + (n * t_k),5) );

  return (n * RELETIVISTIC_F * sqrtA * e_cos_E) / (1.0 - e_cos_E);
}

double Ephemeris::RelTimeRateRate(const double gps_time) const
{
  double A = std::pow(sqrtA,2.0);
  double n_0 = std::sqrt(WGS84_MU / std::pow(A,3.0));
  double t_k = gps_time - t_oe;
  if (t_k > 302400.0) {
    t_k -= 604800.0;
  } else if (t_k < -302400) {
    t_k += 604800.0;
  }  
  double n = n_0 + del_n;
  double E_k = EfromAnomaly(M_0 + (n * t_k),5);

  return ( n * n * RELETIVISTIC_F * std::numbers::e * sqrtA * std::sin(E_k) )
          / std::pow(1.0 - std::numbers::e * std::cos(E_k), 2);
}

void Ephemeris::Randomize()
{
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  M_0 = distribution(internal::random_gen) * (EphemerisUpperLimits.M_0 - EphemerisLowerLimits.M_0) + EphemerisLowerLimits.M_0;
  del_n = distribution(internal::random_gen) * (EphemerisUpperLimits.del_n - EphemerisLowerLimits.del_n) + EphemerisLowerLimits.del_n;
  e = distribution(internal::random_gen) * (EphemerisUpperLimits.e - EphemerisLowerLimits.e) + EphemerisLowerLimits.e;
  sqrtA = distribution(internal::random_gen) * (EphemerisUpperLimits.sqrtA - EphemerisLowerLimits.sqrtA) + EphemerisLowerLimits.sqrtA;
  Omega_0 = distribution(internal::random_gen) * (EphemerisUpperLimits.Omega_0 - EphemerisLowerLimits.Omega_0) + EphemerisLowerLimits.Omega_0;
  i_0 = distribution(internal::random_gen) * (EphemerisUpperLimits.i_0 - EphemerisLowerLimits.i_0) + EphemerisLowerLimits.i_0;
  omega = distribution(internal::random_gen) * (EphemerisUpperLimits.omega - EphemerisLowerLimits.omega) + EphemerisLowerLimits.omega;
  Omega_dot = distribution(internal::random_gen) * (EphemerisUpperLimits.Omega_dot - EphemerisLowerLimits.Omega_dot) + EphemerisLowerLimits.Omega_dot;
  IDOT = distribution(internal::random_gen) * (EphemerisUpperLimits.IDOT - EphemerisLowerLimits.IDOT) + EphemerisLowerLimits.IDOT;
  C_uc = distribution(internal::random_gen) * (EphemerisUpperLimits.C_uc - EphemerisLowerLimits.C_uc) + EphemerisLowerLimits.C_uc;
  C_us = distribution(internal::random_gen) * (EphemerisUpperLimits.C_us - EphemerisLowerLimits.C_us) + EphemerisLowerLimits.C_us;
  C_rc = distribution(internal::random_gen) * (EphemerisUpperLimits.C_rc - EphemerisLowerLimits.C_rc) + EphemerisLowerLimits.C_rc;
  C_rs = distribution(internal::random_gen) * (EphemerisUpperLimits.C_rs - EphemerisLowerLimits.C_rs) + EphemerisLowerLimits.C_rs;
  C_ic = distribution(internal::random_gen) * (EphemerisUpperLimits.C_ic - EphemerisLowerLimits.C_ic) + EphemerisLowerLimits.C_ic;
  C_is = distribution(internal::random_gen) * (EphemerisUpperLimits.C_is - EphemerisLowerLimits.C_is) + EphemerisLowerLimits.C_is;
  t_oe = distribution(internal::random_gen) * (EphemerisUpperLimits.t_oe - EphemerisLowerLimits.t_oe) + EphemerisLowerLimits.t_oe;
}

void Ephemeris::Print() const
{
  std::cout << "M_0:       " << M_0 << '\n';
  std::cout << "del_n:     " << del_n << '\n';
  std::cout << "e:         " << e << '\n';
  std::cout << "sqrtA:     " << sqrtA << '\n';
  std::cout << "Omega_0:   " << Omega_0 << '\n';
  std::cout << "i_0:       " << i_0 << '\n';
  std::cout << "omega:     " << omega << '\n';
  std::cout << "Omega_dot: " << Omega_dot << '\n';
  std::cout << "IDOT:      " << IDOT << '\n';
  std::cout << "C_uc:      " << C_uc << '\n';
  std::cout << "C_us:      " << C_us << '\n';
  std::cout << "C_rc:      " << C_rc << '\n';
  std::cout << "C_rs:      " << C_rs << '\n';
  std::cout << "C_ic:      " << C_ic << '\n';
  std::cout << "C_is:      " << C_is << '\n';
  std::cout << "t_oe:      " << t_oe << '\n';
  std::cout << "IODE:      " << static_cast<int>(IODE) << '\n';
}

} // namespace gps
