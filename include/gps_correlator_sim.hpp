#ifndef SATELLITE_CONSTELLATIONS_INCLUDE_GPS_CORRELATOR_SIM
#define SATELLITE_CONSTELLATIONS_INCLUDE_GPS_CORRELATOR_SIM

#include <complex>
#include <numbers>
#include <random>
#include <iostream>

#include "common_types.hpp"
#include "gps_common.hpp"

namespace Gps {


// Obtained from Scott Martin's dissertation, equation A.3
template<typename FloatType>
void CalcCorrelatorOutput(std::complex<FloatType>& output, const FloatType& chip_error, const FloatType& freq_error,
    const FloatType& phase_error, const FloatType& corr_period, const FloatType& cno)
{
  // Chip Error
  if (std::abs(chip_error) > 1.0) {
    output = 0.0;
    return;
  } else {
    output = (1.0 - std::abs(chip_error));
  }

  // Amplitude
  output *= 2.0 * std::sqrt(cno * corr_period);

  // Frequency Error
  if (freq_error != 0.0) {
    FloatType in_sinc = std::numbers::pi_v<FloatType> * corr_period * freq_error;
    output *= std::sin(in_sinc) / in_sinc;
  }

  // Phase Error
  output *= std::exp(ComplexI<FloatType> * TwoPi<FloatType> * phase_error);
}

template<typename FloatType>
std::complex<FloatType> CalcCorrelatorOutput(const FloatType chip_error, const FloatType freq_error, const FloatType phase_error,
    const FloatType corr_period, const FloatType cno)
{
  std::complex<FloatType> result;
  CalcCorrelatorOutput(result, chip_error, freq_error, phase_error, corr_period, cno);
  return result;
}


// assumes evenly spaced
// template<typename FloatType, BracketContainer Container>
// void CalcCorrelatorOutput(std::complex<FloatType>& output, const Container& chip_errors, const Container& freq_errors,
//     const Container& phase_errors, const FloatType& corr_period, const FloatType& cno)
// {
//   assert(chip_errors.size() == freq_errors.size());
//   assert(freq_errors.size() == phase_errors.size());

//   std::complex<FloatType> section_result;
//   output = 0.;
//   // FloatType section_period = corr_period / static_cast<FloatType>(chip_errors.size());

//   for (std::size_t i = 0; i < chip_errors.size(); i++) {
//     CalcCorrelatorOutput(section_result, chip_errors[i], freq_errors[i], phase_errors[i], corr_period, cno);
//     // CalcCorrelatorOutput(section_result, chip_errors[i], freq_errors[i], phase_errors[i], section_period, cno);
//     output += section_result;
//   }
//   output /= static_cast<FloatType>(chip_errors.size());
// }

template<typename FloatType, BracketContainer Container>
void CalcCorrelatorOutput(std::complex<FloatType>& output, const Container& chip_errors, const Container& freq_errors,
    const Container& phase_errors, const FloatType& corr_period, const FloatType& cno)
{
  assert(chip_errors.size() == freq_errors.size());
  assert(freq_errors.size() == phase_errors.size());

  double chip_average = 0.0;
  double freq_average = 0.0;
  double phase_average = 0.0;
  for (std::size_t i = 0; i < chip_errors.size(); i++) {
    chip_average += chip_errors[i];
    freq_average += freq_errors[i];
    phase_average += phase_errors[i];
    // std::cout << freq_errors[i] << ", ";
  }
  // std::cout << freq_average << std::endl;
  chip_average /= static_cast<FloatType>(chip_errors.size());
  freq_average /= static_cast<FloatType>(chip_errors.size());
  phase_average /= static_cast<FloatType>(chip_errors.size());
  CalcCorrelatorOutput(output, chip_average, freq_average, phase_average, corr_period, cno);
}


template<typename FloatType, BracketContainer Container>
std::complex<FloatType> CalcCorrelatorOutput(const Container& chip_errors, const Container& freq_errors,
    const Container& phase_errors, const FloatType corr_period, const FloatType cno)
{
  std::complex<FloatType> result;
  CalcCorrelatorOutput(result, chip_errors, freq_errors, phase_errors, corr_period, cno);
  return result;
}


template<typename FloatType, bool StoreParams=false>
class CorrelatorSim
{
public:
  CorrelatorSim()
  {}

  ~CorrelatorSim()
  {}

  std::complex<FloatType> Simulate(const FloatType chip_error, const FloatType freq_error, const FloatType phase_error,
      const FloatType corr_period, const FloatType cno) const
  {
    std::complex<FloatType> result;
    CalcCorrelatorOutput(result, chip_error, freq_error, phase_error, corr_period, cno);
    result.real += normal_dist_(internal::random_gen);
    result.imag += normal_dist_(internal::random_gen);
    return result;
  }

  template<BracketContainer Container>
  std::complex<FloatType> Simulate(const Container& chip_errors, const Container& freq_errors,
    const Container& phase_errors, const FloatType corr_period, const FloatType cno) const
  {
    std::complex<FloatType> result;
    CalcCorrelatorOutput<FloatType,Container>(result, chip_errors, freq_errors, phase_errors, corr_period, cno);
    result.real += normal_dist_(internal::random_gen);
    result.imag += normal_dist_(internal::random_gen);
    return result;
  }

  template<bool Exists = StoreParams>
  typename std::enable_if<Exists, std::complex<FloatType>>::type
  Simulate(const FloatType chip_error, const FloatType freq_error, const FloatType phase_error) const
  {
    std::complex<FloatType> result;
    CalcCorrelatorOutput(result, chip_error, freq_error, phase_error, corr_period_, cn_ratio_);
    result.real += normal_dist_(internal::random_gen);
    result.imag += normal_dist_(internal::random_gen);
    return result;
  }

  template<BracketContainer Container, bool Exists = StoreParams>
  typename std::enable_if<Exists, std::complex<FloatType>>::type
  Simulate(const Container& chip_errors, const Container& freq_errors, const Container& phase_errors) const
  {
    std::complex<FloatType> result;
    CalcCorrelatorOutput<FloatType,Container>(result, chip_errors, freq_errors, phase_errors, corr_period_, cn_ratio_);
    result.real += normal_dist_(internal::random_gen);
    result.imag += normal_dist_(internal::random_gen);
    return result;
  }

  template<bool Exists = StoreParams>
  typename std::enable_if<Exists, void>::type SetPeriod(const FloatType corr_period)
  {
    corr_period_ = corr_period;
  }

  template<bool Exists = StoreParams>
  typename std::enable_if<Exists, void>::type SetCNO(const FloatType cn_ratio)
  {
    cn_ratio_ = cn_ratio;
  }

  template<bool Exists = StoreParams>
  typename std::enable_if<Exists, void>::type GetPeriod() const
  {
    return corr_period_;
  }

  template<bool Exists = StoreParams>
  typename std::enable_if<Exists, void>::type GetCNO() const
  {
    return cn_ratio_;
  }

  static void AddNoise(std::complex<FloatType>& output)
  {
    output.real += normal_dist_(internal::random_gen);
    output.imag += normal_dist_(internal::random_gen);
  }

private:
  static std::normal_distribution<FloatType> normal_dist_;

  typename std::enable_if<StoreParams,FloatType>::type corr_period_;
  typename std::enable_if<StoreParams,FloatType>::type cn_ratio_;
};


} // namespace Gps
#endif
