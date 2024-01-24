#ifndef SATELLITE_CONSTELLATIONS_INCLUDE_COMMON_TYPES
#define SATELLITE_CONSTELLATIONS_INCLUDE_COMMON_TYPES

#include <random>
#include <cassert>
#include <complex>
#include <numbers>

// extern std::random_device RANDOM_GENERATOR;

template<class T>
concept BracketContainer = requires (T t, std::size_t index) {
  t.size();
  t[index];
};

template<typename T>
inline constexpr T TwoPi = T(2.) * std::numbers::pi_v<T>;

template<typename T>
inline constexpr std::complex<T> ComplexI = {0.,1.};

template<typename RealType>
void circular_fmod(RealType& x, RealType y)
{
  assert(y > 0.);
  RealType mult = std::floor(x / y);
  x -= mult * y;
}


#endif
