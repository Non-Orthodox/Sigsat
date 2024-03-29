#ifndef SATELLITE_CONSTELLATIONS_INCLUDE_GPS_SIGNAL_GEN
#define SATELLITE_CONSTELLATIONS_INCLUDE_GPS_SIGNAL_GEN

#include <array>
#include <vector>

#include "gps_common.hpp"
#include "gps_lnav_data.hpp"


namespace Gps {
namespace Lnav {

template<typename RealType = double>
struct State
{
  uint8_t subframe = 0; // 0-4
  uint16_t bit = 0; // 0-299
  uint8_t code_cycle = 0; // 0-19
  RealType chip = 0; // [0,1023)
  RealType code_frequency = CA_RATE;
  RealType carrier_frequency = 0.0; // intermediate + doppler
  RealType carrier_phase = 0.0; // radians
};


class SatelliteInfo
{
public:
  SatelliteInfo(const uint8_t prn);
  
  DataFrame& Frame() { return frame_; }
  const std::array<bool,1023>& Code() const { return ca_code_; }
  bool Code(const uint16_t chip_i) const { return ca_code_[chip_i]; }
  
  bool GetMessageBit(const uint8_t subframe_i, const uint16_t bit_i);
  bool Information(const uint8_t subframe_i, const uint16_t bit_i, const uint16_t chip_i);

  void Initialize(uint8_t first_subframe);
private:
  uint8_t GetSubframeIndex(const uint8_t subframe_num);

  uint8_t prn_;

  DataFrame frame_;
  std::array<bool,1023> ca_code_;

  Subframe parity_subframes_ [2];
  uint8_t subframe_nums_ [2];
  //TODO add index for latest subframe? (index of parity_subframes, not subframe index), this will make things more efficient
  //? perhaps make a new class for handling this circular container
};


// function that takes two-buffer set of subframes
//! this function assumes constant carrier frequency and code frequency
template<typename QuantizedType, typename RealType = double>
bool GenSignalWithData(
              uint8_t& subframe,    // 0-4
              uint16_t& bit,        // 0-299
              uint8_t& code_cycle,  // 0-19
              RealType& chip,         // [0,1023)
              RealType& code_frequency,
              RealType& carrier_frequency, // intermediate + doppler
              RealType& carrier_phase, // radians in [0,2*pi]
              SatelliteInfo& sat_info,
              std::complex<QuantizedType>* sample_array, 
              const std::size_t array_size,
              const RealType sample_frequency,
              const QuantizedType amplitude,
              const bool cycle_carryover)
{
  assert(subframe < 5);
  assert(bit < 300);
  assert(code_cycle < 20);

  RealType angular_frequency = carrier_frequency * TwoPi<RealType>;

  RealType prev_chip = cycle_carryover ? 1024.0 : -1.0;
  bool nav_data = sat_info.GetMessageBit(subframe, bit);
	for (uint64_t i = 0; i < array_size; i++) {
    RealType del_t = static_cast<RealType>(i) / sample_frequency;
		RealType current_chip = fmod((del_t * code_frequency) + chip, 1023.0);

    // Update data indices if needed
    if (prev_chip > current_chip) {
      code_cycle++;
      if (code_cycle == 20) {
        code_cycle = 0;
        bit++;
        if (bit == 300) {
          bit = 0;
          subframe++;
          if (subframe == 5) {
            subframe = 0;
          }
        }
        nav_data = sat_info.GetMessageBit(subframe, bit);
      }
    }

		// sample_array[i] = ( (sat_info.Code(static_cast<uint16_t>(current_chip)) ^ nav_data) ? amplitude : -amplitude )
    //                   * std::exp( ComplexI<RealType> * ((angular_frequency * del_t) + carrier_phase) );
    sample_array[i] = amplitude * static_cast<std::complex<QuantizedType>>(
                      ( (sat_info.Code(static_cast<uint16_t>(current_chip)) ^ nav_data) ? 1. : -1. )
                      * std::exp( ComplexI<RealType> * ((angular_frequency * del_t) + carrier_phase) )
                    );
    prev_chip = current_chip;
	}

  RealType del_t = static_cast<RealType>(array_size) / sample_frequency;
	chip = fmod((del_t * code_frequency) + chip, 1023.0);
  carrier_phase = fmod((angular_frequency * del_t) + carrier_phase, TwoPi<RealType>);

  return (prev_chip > chip) ? true : false;
}


template<typename QuantizedType, typename RealType = double>
bool GenSignalWithData(
              State<RealType>& signal_state,
              SatelliteInfo& sat_info,
              std::complex<QuantizedType>* sample_array, 
              const std::size_t array_size,
              const RealType sample_frequency,
              const QuantizedType amplitude,
              const bool cycle_carryover)
{
  return GenSignalWithData<QuantizedType,RealType>(signal_state.subframe,
                    signal_state.bit,
                    signal_state.code_cycle,
                    signal_state.chip,
                    signal_state.code_frequency,
                    signal_state.carrier_frequency,
                    signal_state.carrier_phase,
                    sat_info,
                    sample_array,
                    array_size,
                    sample_frequency,
                    amplitude,
                    cycle_carryover);
}


template<typename QuantizedType, typename RealType = double>
std::vector<bool> GenBasebandSignalsWithData(
              std::vector<State<RealType>>& signal_states,
              std::vector<SatelliteInfo>& sat_info,
              std::complex<QuantizedType>* sample_array, 
              const std::size_t array_size,
              const RealType sample_frequency,
              const RealType amplitude,
              const std::vector<bool> cycle_carryovers)
{
  // std::vector<RealType> angular_frequencies(signal_states.size());
  std::vector<RealType> prev_chips(signal_states.size());
  std::vector<bool> nav_data(signal_states.size());

  for (std::size_t i = 0; i < signal_states.size(); i++) {
    assert(signal_states[i].subframe < 5);
    assert(signal_states[i].bit < 300);
    assert(signal_states[i].code_cycle < 20);

    // angular_frequencies[i] = signal_states[i].carrier_frequency * TwoPi<RealType>;
    prev_chips[i] = cycle_carryovers[i] ? 1024.0 : -1.0;
    nav_data[i] = sat_info[i].GetMessageBit(signal_states[i].subframe, signal_states[i].bit);
  }

  std::vector<RealType> current_chips(signal_states.size());
  for (uint64_t k = 0; k < array_size; k++) {
    RealType del_t = static_cast<RealType>(k) / sample_frequency;

    // For each PRN
    std::complex<RealType> sample = 0.0;
    for (std::size_t i = 0; i < signal_states.size(); i++) {
		  current_chips[i] = fmod((del_t * signal_states[i].code_frequency) + signal_states[i].chip, 1023.0);

      // Update data indices if needed
      if (prev_chips[i] > current_chips[i]) {
        signal_states[i].code_cycle++;
        if (signal_states[i].code_cycle == 20) {
          signal_states[i].code_cycle = 0;
          signal_states[i].bit++;
          if (signal_states[i].bit == 300) {
            signal_states[i].bit = 0;
            signal_states[i].subframe++;
            if (signal_states[i].subframe == 5) {
              signal_states[i].subframe = 0;
            }
          }
          nav_data[i] = sat_info[i].GetMessageBit(signal_states[i].subframe, signal_states[i].bit);
        }
      }

      sample += amplitude * ( (sat_info[i].Code(static_cast<uint16_t>(current_chips[i])) ^ nav_data[i]) ? 1. : -1. );
              // * std::exp( ComplexI<RealType> * ((angular_frequencies[i] * del_t) + signal_states[i].carrier_phase) );
      prev_chips[i] = current_chips[i];
    }
    sample_array[k] = static_cast<std::complex<QuantizedType>>(sample);
  }

  RealType del_t = static_cast<RealType>(array_size) / sample_frequency;
  std::vector<bool> next_carryovers(signal_states.size());

  for (std::size_t i = 0; i < signal_states.size(); i++) {
	  signal_states[i].chip = fmod((del_t * signal_states[i].code_frequency) + signal_states[i].chip, 1023.0);
    next_carryovers[i] = (prev_chips[i] > signal_states[i].chip) ? true : false;
    // signal_states[i].carrier_phase = fmod((angular_frequencies[i] * del_t) + signal_states[i].carrier_phase, TwoPi<RealType>);
  }
  return next_carryovers;
}


} // namespace Lnav
} // namespace Gps

#endif