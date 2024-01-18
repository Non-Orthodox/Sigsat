#ifndef INCLUDE_BINARY_OPS_H
#define INCLUDE_BINARY_OPS_H

#include <cstdint>
#include <array>

bool bitVal(const uint32_t& num, const uint8_t pos);
void bitSet(uint32_t& num, const uint8_t pos);
void bitUnset(uint32_t& num, const uint8_t pos);
void bitToggle(uint32_t& num, const uint8_t pos);
void bitEqu(uint16_t& num, const uint8_t pos, bool val);
void bitEqu(uint32_t& num, const uint8_t pos, bool val);

// prints MSB to LSB
void printBinary(const uint8_t num);
void printBinary(const uint16_t num);
void printBinary(const uint32_t num);

template<int S>
bool multiXOR(const uint32_t& num, const std::array<uint8_t,S> positions)
{
  bool result = bitVal(num, positions[0]);
  for (uint8_t i = 1; i < S; i++) {
    result ^= bitVal(num, positions[i]);
  }
  return result;
}

template<int S>
bool multiXOR(const uint32_t& num, const uint8_t positions[])
{
  bool result = bitVal(num, positions[0]);
  for (uint8_t i = 1; i < S; i++) {
    result ^= bitVal(num, positions[i]);
  }
  return result;
}

#endif