#include <cstdint>
#include <array>
#include <cassert>

#include "binary_ops.hpp"
#include "gps_common.hpp"

namespace Gps
{

namespace internal
{
  std::random_device random_gen;
}


void GenCA(std::array<bool,1023>* const sequence, const uint8_t prn)
{
  assert( !((prn < 1) || (prn > 32)) );

  static constexpr uint8_t g2_out_taps [32][2] = {
    /*01*/{2,6},
    /*02*/{3,7},
    /*03*/{4,8},
    /*04*/{5,9},
    /*05*/{1,9},
    /*06*/{2,10},
    /*07*/{1,8},
    /*08*/{2,9},
    /*09*/{3,10},
    /*10*/{2,3},
    /*11*/{3,4},
    /*12*/{5,6},
    /*13*/{6,7},
    /*14*/{7,8},
    /*15*/{8,9},
    /*16*/{9,10},
    /*17*/{1,4},
    /*18*/{2,5},
    /*19*/{3,6},
    /*20*/{4,7},
    /*21*/{5,8},
    /*22*/{6,9},
    /*23*/{1,3},
    /*24*/{4,6},
    /*25*/{5,7},
    /*26*/{6,8},
    /*27*/{7,9},
    /*28*/{8,10},
    /*29*/{1,6},
    /*30*/{2,7},
    /*31*/{3,8},
    /*32*/{4,9}
  };

  // Linear-feedback shift registers
  uint16_t G1 = 0xFFFF;
  uint16_t G2 = 0xFFFF;

  uint8_t taps1 [2] = {2,9}; // 3,10
  uint8_t taps2 [6] = {1,2,5,7,8,9}; // 2,3,6,8,9,10

  for (std::size_t i = 0; i < 1023; i++) {
    // set value in sequence
    sequence->operator[](i) = bitVal<true>(G1,9) 
                            ^ bitVal<true>(G2,g2_out_taps[prn-1][0]-1)
                            ^ bitVal<true>(G2,g2_out_taps[prn-1][1]-1);

    // shift the registers and set first bits
    bool feedback1 = multiXOR<2,true>(G1,taps1);
    bool feedback2 = multiXOR<6,true>(G2,taps2);
    G1 <<= 1;
    G2 <<= 1;
    bitEqu<true>(G1,0,feedback1);
    bitEqu<true>(G2,0,feedback2);
  }
}


} // namespace gps
