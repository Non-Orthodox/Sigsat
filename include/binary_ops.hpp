#ifndef INCLUDE_BINARY_OPS_H
#define INCLUDE_BINARY_OPS_H

#include <cstdint>
#include <array>
#include <cassert>


template<bool LsbIsZero = true>
bool bitVal(const uint32_t& num, const uint8_t pos)
{
  assert(!(pos > 31));

  if constexpr (LsbIsZero) {
    return (num >> pos) & 1;
  } else {
    return num & (0x80000000 >> pos);
  }
}


template<bool LsbIsZero = true>
void bitSet(uint32_t& num, const uint8_t pos)
{
  assert(!(pos > 31));
  if constexpr (LsbIsZero) {
    num |= 1 << pos;
  } else {
    num |= (0x80000000 >> pos);
  }
}


void bitUnset(uint32_t& num, const uint8_t pos);
void bitToggle(uint32_t& num, const uint8_t pos);


template<bool LsbIsZero = true>
void bitEqu(uint16_t& num, const uint8_t pos, bool val)
{
  assert(!(pos > 15));
  if constexpr (LsbIsZero) {
    num ^= (-(uint16_t)val ^ num) & (1 << pos);
  } else {
    num ^= (-(uint16_t)val ^ num) & (0x8000 >> pos);
  }
}


template<bool LsbIsZero = true>
void bitEqu(uint32_t& num, const uint8_t pos, bool val)
{
  assert(!(pos > 31));
  if constexpr (LsbIsZero) {
    num ^= (-(uint32_t)val ^ num) & (1 << pos);
  } else {
    num ^= (-(uint32_t)val ^ num) & (0x80000000 >> pos);
  }
}


// prints MSB to LSB
void printBinary(const uint8_t num);
void printBinary(const uint16_t num);
void printBinary(const uint32_t num);

template<int Size, bool IsLsbZero = true>
bool multiXOR(const uint32_t& num, const std::array<uint8_t,Size> positions)
{
  bool result = bitVal<IsLsbZero>(num, positions[0]);
  for (uint8_t i = 1; i < Size; i++) {
    result ^= bitVal<IsLsbZero>(num, positions[i]);
  }
  return result;
}

template<int Size, bool IsLsbZero = true>
bool multiXOR(const uint32_t& num, const uint8_t positions[])
{
  bool result = bitVal<IsLsbZero>(num, positions[0]);
  for (uint8_t i = 1; i < Size; i++) {
    result ^= bitVal<IsLsbZero>(num, positions[i]);
  }
  return result;
}

#endif