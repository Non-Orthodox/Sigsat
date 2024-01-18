#include <cstdint>
#include <iostream>
#include <cassert>

#include "binary_ops.hpp"

// static_assert(sizeof(uint8_t) == 1);
// static_assert(sizeof(uint16_t) == 2);
// static_assert(sizeof(uint32_t) == 4);
// static_assert(sizeof(uint64_t) == 8);
// static_assert(sizeof(int8_t) == 1);
// static_assert(sizeof(int16_t) == 2);
// static_assert(sizeof(int32_t) == 4);
// static_assert(sizeof(int64_t) == 8);

template<typename IntType>
bool bitVal(const IntType& num, const uint8_t pos)
{
  assert(pos > (sizeof(IntType)*8-1));
  return (num >> pos) & 1;
}

bool bitVal(const uint32_t& num, const uint8_t pos)
{
  assert(!(pos > 31));
  return (num >> pos) & 1;
}

template<typename IntType>
void bitSet(IntType& num, const uint8_t pos)
{
  assert(pos > (sizeof(IntType)*8-1));
  num |= 1 << pos;
}

void bitSet(uint32_t& num, const uint8_t pos)
{
  assert(!(pos > 31));
  num |= 1 << pos;
}

template<typename IntType>
void bitUnset(IntType& num, const uint8_t pos)
{
  assert(pos > (sizeof(IntType)*8-1));
  num &= ~(1 << pos);
}

void bitUnset(uint32_t& num, const uint8_t pos)
{
  assert(!(pos > 31));
  num &= ~(1 << pos);
}

template<typename IntType>
void bitToggle(IntType& num, const uint8_t pos)
{
  assert(pos > (sizeof(IntType)*8-1));
  num ^= 1 << pos;
}

void bitToggle(uint32_t& num, const uint8_t pos)
{
  assert(!(pos > 31));
  num ^= 1 << pos;
}

void bitEqu(uint16_t& num, const uint8_t pos, bool val)
{
  assert(!(pos > 15));
  num ^= (-(uint16_t)val ^ num) & (1 << pos);
}

void bitEqu(uint32_t& num, const uint8_t pos, bool val)
{
  assert(!(pos > 31));
  num ^= (-(uint32_t)val ^ num) & (1 << pos);
}

void printBinary(const uint8_t num)
{
  for (uint8_t i = 0; i < 8; i++) {
    std::cout << bitVal(num, 7 - i);
  }
  std::cout << std::endl;
}

void printBinary(const uint16_t num)
{
  for (uint8_t i = 0; i < 16; i++) {
    std::cout << bitVal(num, 15 - i);
  }
  std::cout << std::endl;
}

void printBinary(const uint32_t num)
{
  for (uint8_t i = 0; i < 32; i++) {
    std::cout << bitVal(num, 31 - i);
  }
  std::cout << std::endl;
}