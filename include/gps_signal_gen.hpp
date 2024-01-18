#ifndef SATELLITE_CONSTELLATIONS_INCLUDE_GPS_SIGNAL_GEN
#define SATELLITE_CONSTELLATIONS_INCLUDE_GPS_SIGNAL_GEN

#include <array>

#include "gps_lnav_data.hpp"


namespace Gps {
namespace Lnav {

template<typename ScalarType = double>
struct LnavState
{
  uint8_t subframe; // 0-4
  uint16_t bit; // 0-299
  uint8_t code_cycle; // 0-19
  ScalarType chip; // [0,1023)
  ScalarType code_frequency;
  ScalarType carrier_frequency; // intermediate + doppler
  ScalarType carrier_phase; // radians
};


class SatelliteLnavInfo
{
public:
  SatelliteLnavInfo(const uint8_t prn);
  
  DataFrame& Frame() { return frame_; }
  const std::array<bool,1023>& Code() const { return ca_code_; }
  bool Code(const uint16_t chip_i) const { return ca_code_[chip_i]; }
  
  bool GetMessageBit(const uint8_t subframe_i, const uint8_t bit_i);
  bool Information(const uint8_t subframe_i, const uint8_t bit_i, const uint16_t chip_i);

private:
  uint8_t GetSubframeIndex(const uint8_t subframe_num);
  void Initialize(uint8_t first_subframe);

  uint8_t prn_;

  DataFrame frame_;
  std::array<bool,1023> ca_code_;

  Subframe parity_subframes_ [2];
  uint8_t subframe_nums_ [2];
  //TODO maybe add index for latest subframe?
  //? perhaps make a new class for handling this circular container
};


// function that takes two-buffer set of subframes
template<typename QuantizedType, typename RealType = double>
void GenSignalWithData(
              uint8_t& subframe,    // 0-4
              uint16_t& bit,        // 0-299
              uint8_t& code_cycle,  // 0-19
              RealType& chip,         // [0,1023)
              RealType& code_frequency,
              RealType& carrier_frequency, // intermediate + doppler
              RealType& carrier_phase, // radians
              SatelliteLnavInfo& sat_info,
              std::complex<QuantizedType>* const sample_array, 
              const std::size_t array_size,
              const RealType sample_frequency,
              const QuantizedType amplitude)
{
  assert(subframe <= 5);
  assert(bit <= 300);
  assert(code_cycle <= 20);

  RealType angular_frequency = carrier_frequency * TwoPi<RealType>;
  RealType prev_chip = -1.0;
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

		sample_array[i] = ( (sat_info.Code(static_cast<uint16_t>(current_chip)) ^ nav_data) ? amplitude : -amplitude )
                      * std::exp( ComplexI<RealType> * ((angular_frequency * del_t) + carrier_phase) );
    prev_chip = current_chip;
	}
  RealType del_t = static_cast<RealType>(array_size) / sample_frequency;
	chip = fmod((del_t * code_frequency) + chip, 1023.0);
  carrier_phase = fmod((angular_frequency * del_t) + carrier_phase, TwoPi<RealType>);
}


template<typename ScalarType>
void GenSignalWithData(
              LnavState<ScalarType>& signal_state,
              const SatelliteLnavInfo& sat_info,
              std::complex<ScalarType>* const sample_array, 
              const std::size_t array_size,
              const ScalarType sample_frequency,
              const ScalarType amplitude)
{
  GenSignalWithData(signal_state.subframe,
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
                    amplitude);
}


} // namespace Lnav
} // namespace Gps

#endif