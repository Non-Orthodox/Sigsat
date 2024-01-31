#include <array>
#include <algorithm>
#include <iostream>

#include "gps_lnav_data.hpp"


namespace Gps {
namespace Lnav {


void TLM(Word& word, const uint16_t tlm_message, const bool integrity_status_flag)
{
  word.Set(0);
  word.Set(4);
  word.Set(6);
  word.Set(7);
  word.SegmentSet(8, tlm_message, 0, 13);
  word.Set(22, integrity_status_flag);
}

Word TLM(const uint16_t tlm_message, const bool integrity_status_flag)
{
  Word word;
  TLM(word, tlm_message, integrity_status_flag);
  return word;
}

void HOW(Word& word, const uint32_t full_tow, const bool alert_flag,
        const bool anti_spoof_flag, const uint8_t subframe_id)
{
  assert( !((subframe_id < 1) || (subframe_id > 5)) );
  // std::cerr << "gps_common::HOW Error: subframe_num " << subframe_id << " invalid." << std::endl;
  
  word.SegmentSet(0, full_tow, 2, 18);
  word.Set(17, alert_flag);
  word.Set(18, anti_spoof_flag);
  word.SegmentSet(19, subframe_id, 0, 2);
}

Word HOW(const uint32_t full_tow, const bool alert_flag,
        const bool anti_spoof_flag, const uint8_t subframe_id)
{
  Word word;
  HOW(word, full_tow, alert_flag, anti_spoof_flag, subframe_id);
  return word;
}


//---------------------------------- Word ----------------------------------
void Word::SegmentSet(uint8_t loc, const uint32_t& val, uint8_t start, uint8_t end)
{
  // TODO assertions
  for (uint8_t i = 0; i <= end-start; i++) {
    Set(loc+i, val, end-i);
  }
}

void Word::SegmentSet(uint8_t loc, const uint32_t& val, uint8_t end)
{
  for (uint8_t i = 0; i <= end; i++) {
    Set(loc+i, val, end-i);
  }
}

Word Word::Parity(bool D29, bool D30) const
{
  uint32_t word = bits_;
  bitEqu(word, 24, (D29 ^ multiXOR<14>(word, parity_array_25)));
  bitEqu(word, 25, (D30 ^ multiXOR<14>(word, parity_array_26)));
  bitEqu(word, 26, (D29 ^ multiXOR<14>(word, parity_array_27)));
  bitEqu(word, 27, (D30 ^ multiXOR<14>(word, parity_array_28)));
  bitEqu(word, 28, (D30 ^ multiXOR<15>(word, parity_array_29)));
  bitEqu(word, 29, (D29 ^ multiXOR<13>(word, parity_array_30)));
  for (uint8_t i = 0; i < 24; i++) {
    bitEqu(word, i, bitVal(word,i) ^ D30);
  }
  return word;
}

void Word::Print() const
{
  for (uint8_t i = 0; i < 30; i++) {
    std::cout << bitVal(bits_, i);
  }
  std::cout << '\n';
}


//---------------------------------- Subframe ----------------------------------
Word& Subframe::operator[](const uint8_t index)
{
  assert( !(index > 9) );
  return words_[index];
}

Word Subframe::operator[](const uint8_t index) const
{
  assert( !(index > 9) );
  return words_[index];
}

bool Subframe::Bit(const uint16_t index) const
{
  return words_[index / 30].Bit(index % 30); 
}

void Subframe::Print() const
{
  for (uint8_t i = 0; i < 10; i++) {
    words_[i].Print();
  }
}


//---------------------------------- DataFrame ----------------------------------
Subframe& DataFrame::operator[](uint8_t index)
{
  return subframes_[index];
}

Subframe DataFrame::operator[](uint8_t index) const
{
  return subframes_[index];
}

// should be done every new frame
void DataFrame::TimeIncrement()
{
  tow_ += 20;
  if (tow_ > 403199) {
    tow_ = tow_ % 403200;
    week_++;
  }
}

Subframe DataFrame::ParityFrame(uint8_t sf)
{
  Subframe result = subframes_[sf];

  //                           1,3,5,6,7,9,10,14,15,16,17,18,21,22
  static uint8_t arr29 [14] = {0,2,4,5,6,8, 9,13,14,15,16,17,20,21};
  //                           3,5,6,8,9,10,11,13,15,19,22,24
  static uint8_t arr30 [12] = {2,4,5,7,8, 9,10,12,14,18,21,23};

  for (uint8_t w = 0; w < 10; w++) {
    // setting bearing bits
    if ((w == 1) || (w == 9)) {
      result[w].Set(23, (D30_ ^ multiXOR<14>( result[w].Val(), arr29 )));
      result[w].Set(22, (D29_ ^ multiXOR<12>( result[w].Val(), arr30 )));
    }
    result[w] = result[w].Parity(D29_,D30_);
    D29_ = result[w].Bit(28);
    D30_ = result[w].Bit(29);
  }
  return result;
}

void DataFrame::SetSubframe(uint8_t sf) // 0-indexed
{
  switch(sf) {
    case 0:
      SetSubframe1();
      return;
    case 1:
      SetSubframe2();
      return;
    case 2:
      SetSubframe3();
      return;
    case 3:
      SetSubframe4();
      return;
    case 4:
      SetSubframe5();
      return;
  }
}

void DataFrame::SetSubframes()
{
  SetSubframe1();
  SetSubframe2();
  SetSubframe3();
  SetSubframe4();
  SetSubframe5();
}

void DataFrame::Preamble(const uint8_t sf_i)
{
  TLM(subframes_[sf_i][0], 0, integrity_status_flag_);
  HOW(subframes_[sf_i][1], tow_ + (sf_i * 4), alert_flag_, anti_spoof_flag_, sf_i + 1);
}

void DataFrame::SetSubframe1()
{
  Preamble(0);

  subframes_[0][2].SegmentSet(0,  week_, 0, 9);
  subframes_[0][2].SegmentSet(10, l2_flag_, 0, 1);
  subframes_[0][2].SegmentSet(12, URA_, 0, 3);
  subframes_[0][2].SegmentSet(16, health_, 0, 5);
  subframes_[0][2].SegmentSet(22, clock_data_.IODC, 8, 9);
  
  // Words 4,5,6 are reserved

  subframes_[0][6].SegmentSet(16, T_GD(), 0, 7);

  subframes_[0][7].SegmentSet(0, clock_data_.IODC, 0, 7);
  subframes_[0][7].SegmentSet(8, t_oc(), 0, 15);

  subframes_[0][8].SegmentSet(0, a_f2(), 0, 7);
  subframes_[0][8].SegmentSet(8, a_f1(), 0, 15);

  subframes_[0][9].SegmentSet(0, a_f0(), 0, 21);
}

void DataFrame::SetSubframe2()
{
  Preamble(1);

  subframes_[1][2].SegmentSet(0, ephemeris_.IODE, 0, 7);
  subframes_[1][2].SegmentSet(8, ephemeris_.C_rs, 0, 15);

  subframes_[1][3].SegmentSet(0, del_n(), 0, 15);
  int32_t M = M_0();
  subframes_[1][3].SegmentSet(16, M, 24, 31);

  subframes_[1][4].SegmentSet(0, M, 0, 23);

  subframes_[1][5].SegmentSet(0, C_uc(), 0, 15);
  int32_t e_bin = e();
  subframes_[1][5].SegmentSet(16, e_bin, 24, 31);

  subframes_[1][6].SegmentSet(0, e_bin, 0, 23);

  subframes_[1][7].SegmentSet(0, C_us(), 0, 15);
  int32_t sqrtA_bin = sqrtA();
  subframes_[1][7].SegmentSet(16, sqrtA_bin, 24, 31);

  subframes_[1][8].SegmentSet(0, sqrtA_bin, 0, 23);

  subframes_[1][9].SegmentSet(0, t_oe(), 0, 15);
  subframes_[1][9].Set(16, fit_interval_flag_);
  subframes_[1][9].SegmentSet(17, AODO_, 0, 4);
}

//TODO check
void DataFrame::SetSubframe3()
{
  Preamble(2);

  subframes_[2][2].SegmentSet(0, C_ic(), 0, 15);
  int32_t Omega0_bin = Omega_0();
  subframes_[2][2].SegmentSet(16, Omega0_bin, 24, 31);

  subframes_[2][3].SegmentSet(0, Omega0_bin, 0, 23);

  subframes_[2][4].SegmentSet(0, C_is(), 0, 15);
  int32_t i0_bin = i_0();
  subframes_[2][4].SegmentSet(16, i0_bin, 24, 31);

  subframes_[2][5].SegmentSet(0, i0_bin, 0, 23);

  subframes_[2][6].SegmentSet(0, C_rc(), 0, 15);
  int32_t omega_bin = omega();
  subframes_[2][6].SegmentSet(16, omega_bin, 24, 31);
  
  subframes_[2][7].SegmentSet(0, omega_bin, 0, 23);

  subframes_[2][8].SegmentSet(0, Omega_dot(), 0, 23);

  subframes_[2][9].SegmentSet(0, ephemeris_.IODE, 0, 7);
  subframes_[2][9].SegmentSet(8, IDOT(), 0, 13);
}

void DataFrame::SetSubframe4()
{
  constexpr std::array<uint8_t,13> reserved_pages = {1,6,11,12,14,15,16,19,20,21,22,23,24};
  constexpr std::array<uint8_t,8> almanac_pages = {2,3,4,5,7,8,9,10};
  constexpr std::array<uint8_t,25> subframe4_ids = 
    {57,25,26,27,28,57,29,30,31,32,57,62,52,53,54,57,55,56,58,59,57,60,61,62,63};

  Preamble(3);

  // First 2 bits data ID, 10 is only valid value
  subframes_[3][2].Set(0);
  subframes_[3][2].SegmentSet(2,subframe4_ids[page_],0,5); // TODO cast id to uint32_t?
  
  if (std::find(reserved_pages.begin(), reserved_pages.end(), page_+1) != reserved_pages.end()) {
    return;
  }

  if (std::find(almanac_pages.begin(), almanac_pages.end(), page_+1) != almanac_pages.end()) {
    // write almanac data
    return;
  }

  if ((page_+1) == 13) {
    // NMCT (nav message correction table)
    return;
  }

  if ((page_+1) == 17) {
    // "special messages"
    return;
  }

  if ((page_+1) == 18) {
    // ionospheric and UTC data
    return;
  }

  if ((page_+1) == 25) {
    // A-S flags, SV configs, health
    return;
  }
}

void DataFrame::SetSubframe5()
{
  constexpr std::array<uint8_t,25> subframe5_ids = 
  {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,51};

  Preamble(4);

  subframes_[3][2].Set(0);
  subframes_[3][2].SegmentSet(2,subframe5_ids[page_],0,5); // TODO cast id to uint32_t?

  if ((page_ >= 0) || (page_ < 24)) {
    // almanac data
    return;
  }

  if (page_ == 24) {
    // sv health
    return;
  }

}

int8_t DataFrame::T_GD()
{
  return ParamToBinary(clock_data_.T_GD, ClockDataScaleFactors.T_GD);
}

uint16_t DataFrame::t_oc()
{
  return ParamToBinary(clock_data_.t_oc, ClockDataScaleFactors.t_oc);
}

int32_t DataFrame::a_f0()
{
  return ParamToBinary(clock_data_.a_f0, ClockDataScaleFactors.a_f0);
}

int16_t DataFrame::a_f1()
{
  return ParamToBinary(clock_data_.a_f1, ClockDataScaleFactors.a_f1);
}

int8_t DataFrame::a_f2()
{
  return ParamToBinary(clock_data_.a_f2, ClockDataScaleFactors.a_f2);
}

int16_t DataFrame::C_rs()
{
  return ParamToBinary(ephemeris_.C_rs, EphemerisScaleFactors.C_rs);
}

int16_t DataFrame::del_n()
{
  return ParamToBinary(ephemeris_.del_n, EphemerisScaleFactors.del_n);
}

int32_t DataFrame::M_0()
{
  return ParamToBinary(ephemeris_.M_0, EphemerisScaleFactors.M_0);
}

int16_t DataFrame::C_uc()
{
  return ParamToBinary(ephemeris_.C_uc, EphemerisScaleFactors.C_uc);
}

int32_t DataFrame::e()
{
  return ParamToBinary(ephemeris_.e, EphemerisScaleFactors.e);
}

int16_t DataFrame::C_us()
{
  return ParamToBinary(ephemeris_.C_us, EphemerisScaleFactors.C_us);
}

uint32_t DataFrame::sqrtA()
{
  return ParamToBinary(ephemeris_.sqrtA, EphemerisScaleFactors.sqrtA);
}

uint16_t DataFrame::t_oe()
{
  return ParamToBinary(ephemeris_.t_oe, EphemerisScaleFactors.t_oe);
}

int16_t DataFrame::C_ic()
{
  return ParamToBinary(ephemeris_.C_ic, EphemerisScaleFactors.C_ic);
}

int32_t DataFrame::Omega_0()
{
  return ParamToBinary(ephemeris_.Omega_0, EphemerisScaleFactors.Omega_0);
}

int16_t DataFrame::C_is()
{
  return ParamToBinary(ephemeris_.C_is, EphemerisScaleFactors.C_is);
}

int32_t DataFrame::i_0()
{
  return ParamToBinary(ephemeris_.i_0, EphemerisScaleFactors.i_0);
}

int16_t DataFrame::C_rc()
{
  return ParamToBinary(ephemeris_.C_rc, EphemerisScaleFactors.C_rc);
}

int32_t DataFrame::omega()
{
  return ParamToBinary(ephemeris_.omega, EphemerisScaleFactors.omega);
}

int32_t DataFrame::Omega_dot()
{
  return ParamToBinary(ephemeris_.Omega_dot, EphemerisScaleFactors.Omega_dot);
}

int16_t DataFrame::IDOT()
{
  return ParamToBinary(ephemeris_.IDOT, EphemerisScaleFactors.IDOT);
}

void DataFrame::RandomizeParams()
{
  clock_data_.Randomize();
  ephemeris_.Randomize();
  clock_data_.t_oc = ephemeris_.t_oe;

  ephemeris_.IODE = 241;
  clock_data_.IODC = ephemeris_.IODE;
}


} // namespace Lnav
} // namespace Gps