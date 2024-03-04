#ifndef SATELLITE_CONSTELLATIONS_INCLUDE_GPS_COMMON
#define SATELLITE_CONSTELLATIONS_INCLUDE_GPS_COMMON

#include <cstdint>
#include <complex>
#include <array>
#include <cmath>
#include <cassert>
#include <ctime>

#include "common_types.hpp"

namespace Gps
{

namespace internal
{
  extern std::random_device random_gen;
}

// These are the exact precision as specified in the IS-GPS-200 documentation
const double LIGHT_SPEED = 299792458.0;
const double PI = 3.1415926535898; // used in ephemeris calculations

const double PI_2 = 2.0 * PI;

const double L1_FREQUENCY = 1.57542e9;
const double L1_ANGULAR_FREQUENCY = TwoPi<double> * L1_FREQUENCY;

const double L2_FREQUENCY = 1.22760e9;
const double L5_FREQUENCY = 1.17645e9;
const double CA_RATE = 1.023e6;
const double P_RATE = 10.23e6;
const double M_RATE = 5.115e6;
const double L5_CODE_RATE = 10.23e6;

const double DATA_BIT_RATE = 50.0;
const double DATA_BIT_PERIOD = 0.02;

const double CHIPS_PER_METER = CA_RATE / LIGHT_SPEED;
const double L1_RADIANS_PER_METER = L1_ANGULAR_FREQUENCY / LIGHT_SPEED;


void GenCA(std::array<bool,1023>* const sequence, const uint8_t prn);


//--------------------------- Correlation ---------------------------
// Container objects must have bracket operator [] for accessing elements
template<template<class> class Container, typename ScalarType>
ScalarType Correlate(const Container<ScalarType>& vec1, const Container<ScalarType>& vec2)
{
  assert(vec1.size() == vec2.size());
  ScalarType result = 0.0;
  for (std::size_t i = 0; i < vec1.size(); i++) {
    result += vec1[i] * vec2[i];
  }
  return result;
}

template<template<class> class Container, typename ScalarType>
std::complex<ScalarType> ComplexCorrelate(const Container<std::complex<ScalarType>>& vec1,
  const Container<ScalarType>& vec2)
{
  assert(vec1.size() == vec2.size());
  std::complex<ScalarType> result = 0.0;
  for (std::size_t i = 0; i < vec1.size(); i++) {
    result += vec1[i] * vec2[i];
  }
  return result;
}

template<template<class> class Container, typename ScalarType>
std::complex<ScalarType> ComplexCorrelate(const Container<ScalarType>& vec1,
  const Container<std::complex<ScalarType>>& vec2)
{
  assert(vec1.size() == vec2.size());
  std::complex<ScalarType> result = 0.0;
  for (std::size_t i = 0; i < vec1.size(); i++) {
    result += vec1[i] * vec2[i];
  }
  return result;
}

template<template<class> class Container, typename ScalarType>
std::complex<double> ComplexCorrelate(const Container<std::complex<ScalarType>>& vec1,
  const Container<std::complex<ScalarType>>& vec2)
{
  assert(vec1.size() == vec2.size());
  std::complex<double> result = 0.0;
  for (std::size_t i = 0; i < vec1.size(); i++) {
    result += vec1[i] * std::conj(vec2[i]);
  }
  return result;
}


//--------------------------- CA Code Sampling ---------------------------
template<typename QuantizedType, typename RealType = double>
void SampleBasebandCa(const std::array<bool,1023>& ca_code, QuantizedType* const sample_array,
              const std::size_t array_size, const RealType sample_frequency, RealType& start_chip,
              const QuantizedType amplitude)
{
	for (uint64_t i = 0; i < array_size; i++) {
		RealType current_chip = fmod((static_cast<RealType>(i) * CA_RATE / sample_frequency) + start_chip, 1023.0);
		sample_array[i] = ca_code[ static_cast<uint16_t>(current_chip) ] ? amplitude : -amplitude;
	}
	start_chip = fmod((static_cast<RealType>(array_size) * CA_RATE / sample_frequency) + start_chip, 1023.0);
}

template<typename QuantizedType, typename RealType = double>
void SampleBasebandCa(const std::array<bool,1023>& ca_code, std::complex<QuantizedType>* const sample_array,
              const std::size_t array_size, const RealType sample_frequency, RealType& start_chip,
              const QuantizedType amplitude)
{
	for (uint64_t i = 0; i < array_size; i++) {
		RealType current_chip = fmod((static_cast<RealType>(i) * CA_RATE / sample_frequency) + start_chip, 1023.0);
		sample_array[i] = ca_code[ static_cast<uint16_t>(current_chip) ] ? amplitude : -amplitude;
	}
	start_chip = fmod((static_cast<RealType>(array_size) * CA_RATE / sample_frequency) + start_chip, 1023.0);
}


// Allows choice of code frequency
template<typename QuantizedType, typename RealType = double>
void SampleBasebandCa(const std::array<bool,1023>& ca_code, QuantizedType* const sample_array,
              const std::size_t array_size, const RealType sample_frequency, const RealType code_frequency,
              RealType& start_chip, const QuantizedType amplitude)
{
	for (uint64_t i = 0; i < array_size; i++) {
		RealType current_chip = fmod((static_cast<RealType>(i) * code_frequency / sample_frequency) + start_chip, 1023.0);
		sample_array[i] = ca_code[ static_cast<uint16_t>(current_chip) ] ? amplitude : -amplitude;
	}
	start_chip = fmod((static_cast<RealType>(array_size) * code_frequency / sample_frequency) + start_chip, 1023.0);
}

template<typename QuantizedType, typename RealType = double>
void SampleBasebandCa(const std::array<bool,1023>& ca_code, std::complex<QuantizedType>* const sample_array,
              const std::size_t array_size, const RealType sample_frequency, const RealType code_frequency,
              RealType& start_chip, const QuantizedType amplitude)
{
	for (uint64_t i = 0; i < array_size; i++) {
		RealType current_chip = fmod((static_cast<RealType>(i) * code_frequency / sample_frequency) + start_chip, 1023.0);
		sample_array[i] = ca_code[ static_cast<uint16_t>(current_chip) ] ? amplitude : -amplitude;
	}
	start_chip = fmod((static_cast<RealType>(array_size) * code_frequency / sample_frequency) + start_chip, 1023.0);
}


// Sample CA code with specified code frequency and intermediate frequency
template<typename QuantizedType, typename RealType = double>
void SampleCa(const std::array<bool,1023>& ca_code, QuantizedType* const sample_array,
              const std::size_t array_size, const RealType sample_frequency,
              const RealType code_frequency, RealType& start_chip, const RealType intermediate_frequency,
              RealType& carrier_phase, const QuantizedType amplitude)
{
  RealType angular_frequency = intermediate_frequency * TwoPi<RealType>;
	for (uint64_t i = 0; i < array_size; i++) {
    RealType del_t = static_cast<RealType>(i) / sample_frequency;
		RealType current_chip = fmod((del_t * code_frequency) + start_chip, 1023.0);
		sample_array[i] = ca_code[ static_cast<uint16_t>(current_chip) ] ? amplitude : -amplitude;
    sample_array[i] *= std::cos((angular_frequency * del_t) + carrier_phase);
	}
  RealType del_t = static_cast<RealType>(array_size) / sample_frequency;
	start_chip = fmod((del_t * code_frequency) + start_chip, 1023.0);
  carrier_phase = fmod((angular_frequency * del_t) + carrier_phase, TwoPi<RealType>);
}

// Note: RealType must be compatible with the complex exponential function
template<typename QuantizedType, typename RealType = double>
void SampleCa(const std::array<bool,1023>& ca_code, std::complex<QuantizedType>* const sample_array, 
              const std::size_t array_size, const RealType sample_frequency,
              const RealType code_frequency, RealType& start_chip, const RealType intermediate_frequency,
              RealType& carrier_phase, const QuantizedType amplitude)
{
  RealType angular_frequency = intermediate_frequency * TwoPi<RealType>;
	for (uint64_t i = 0; i < array_size; i++) {
    RealType del_t = static_cast<RealType>(i) / sample_frequency;
		RealType current_chip = fmod((del_t * code_frequency) + start_chip, 1023.0);
		sample_array[i] = ( ca_code[ static_cast<uint16_t>(current_chip) ] ? amplitude : -amplitude )
                       * std::exp( ComplexI<RealType> * ((angular_frequency * del_t) + carrier_phase) ); 
	}
  RealType del_t = static_cast<RealType>(array_size) / sample_frequency;
	start_chip = fmod((del_t * code_frequency) + start_chip, 1023.0);
  carrier_phase = fmod((angular_frequency * del_t) + carrier_phase, TwoPi<RealType>);
}


} // namespace gps
#endif
