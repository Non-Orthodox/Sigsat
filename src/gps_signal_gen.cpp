
#include "gps_common.hpp"
#include "gps_signal_gen.hpp"

namespace Gps {
namespace Lnav {

SatelliteLnavInfo::SatelliteLnavInfo(const uint8_t prn) : prn_{prn}
{
  // TODO modify data frame?
  Initialize(0);
}


void SatelliteLnavInfo::Initialize(uint8_t first_subframe)
{
  GenCA(&ca_code_,prn_);
  subframe_nums_[0] = first_subframe;
  subframe_nums_[1] = (first_subframe + 1) % 5;
}


uint8_t SatelliteLnavInfo::GetSubframeIndex(const uint8_t subframe_num)
{
  assert(subframe_num < 5);
  if (subframe_nums_[0] == subframe_num) {
    return 0;
  } 
  else if (subframe_nums_[1] == subframe_num) {
    return 1;
  }
  else {
    // find which one to swap
    if (subframe_nums_[1] == ((subframe_nums_[0] + 1) % 5)) {
      subframe_nums_[0] = (subframe_nums_[1] + 1) % 5;
      parity_subframes_[0] = frame_.ParityFrame(subframe_nums_[0]);
      assert(subframe_nums_[0] == subframe_num);
      return 0;
    }
    else {
      subframe_nums_[1] = (subframe_nums_[0] + 1) % 5;
      parity_subframes_[1] = frame_.ParityFrame(subframe_nums_[1]);
      assert(subframe_nums_[1] == subframe_num);
      return 1;
    }
  }
}


bool SatelliteLnavInfo::GetMessageBit(const uint8_t subframe_i, const uint8_t bit_i)
{
  uint8_t subframe_index = GetSubframeIndex(subframe_i);
  return parity_subframes_[subframe_index].Bit(bit_i);
}


bool SatelliteLnavInfo::Information(const uint8_t subframe_i, const uint8_t bit_i, const uint16_t chip_i)
{
  return GetMessageBit(subframe_i, bit_i) ^ ca_code_[chip_i];
}


} // namespace Lnav
} // namespace Gps