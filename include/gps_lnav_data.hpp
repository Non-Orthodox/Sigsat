#ifndef SATELLITE_CONSTELLATIONS_INCLUDE_GPS_LNAV_DATA
#define SATELLITE_CONSTELLATIONS_INCLUDE_GPS_LNAV_DATA

#include "binary_ops.hpp"
#include "gps_ephemeris.hpp"

namespace Gps {
namespace Lnav {

//                                       {1,2,3,5,6,10,11,12,13,14,17,18,20,23}
constexpr uint8_t parity_array_25 [14] = {0,1,2,4,5, 9,10,11,12,13,16,17,19,22};
//                                       {2,3,4,6,7,11,12,13,14,15,18,19,21,24}
constexpr uint8_t parity_array_26 [14] = {1,2,3,5,6,10,11,12,13,14,17,18,20,23};
//                                       {1,3,4,5,7, 8,12,13,14,15,16,19,20,22}
constexpr uint8_t parity_array_27 [14] = {0,2,3,4,6, 7,11,12,13,14,15,18,19,21};
//                                       {2,4,5,6,8, 9,13,14,15,16,17,20,21,23}
constexpr uint8_t parity_array_28 [14] = {1,3,4,5,7, 8,12,13,14,15,16,19,20,22};
//                                       {1,3,5,6,7, 9,10,14,15,16,17,18,21,22,24}
constexpr uint8_t parity_array_29 [15] = {0,2,4,5,6, 8, 9,13,14,15,16,17,20,21,23};
//                                       {3,5,6,8,9,10,11,13,15,19,22,23,24}
constexpr uint8_t parity_array_30 [13] = {2,4,5,7,8, 9,10,12,14,18,21,22,23};


class Word
{
private:
  uint32_t bits_ {0};

public:
  Word() {}
  Word(uint32_t bits) : bits_{bits} {}
  ~Word() {}
  
  // Convenient Operators
  Word& operator=(const Word other) { this->bits_ = other.bits_; return *this; }
  Word& operator=(const uint32_t bits) { bits_ = bits; return *this; }
  
  // Reading Data (cannot write with these)
  bool operator[](const uint8_t pos) const { return bitVal(bits_, pos); }
  bool Bit(const uint8_t pos) const { return bitVal(bits_, pos); }

  // Changing Bits
  void Set(uint8_t pos) { bitSet(bits_, pos); }
  void Set(uint8_t pos, bool val) { bitEqu(bits_,pos, val); }
  void Set(uint8_t pos, uint32_t val, uint8_t valpos) { bitEqu(bits_, pos, bitVal(val, valpos)); }
  void SetVal(uint32_t val) { bits_=val; }
  void SegmentSet(uint8_t loc, const uint32_t& val, uint8_t start, uint8_t end);
  void SegmentSet(uint8_t loc, const uint32_t& val, uint8_t end);

  // void Parity(bool D29, bool D30);
  Word Parity(bool D29, bool D30) const;

  uint32_t Val() const { return bits_; }
  void Reset() { bits_ = 0; }
  void Print() const; // prints MSB to LSB of message
};


void TLM(Word& word, const uint16_t tlm_message, const bool integrity_status_flag);

Word TLM(const uint16_t tlm_message, const bool integrity_status_flag);

void HOW(Word& word, const uint32_t full_tow, const bool alert_flag,
        const bool anti_spoof_flag, const uint8_t subframe_id);

Word HOW(const uint32_t full_tow, const bool alert_flag,
        const bool anti_spoof_flag, const uint8_t subframe_id);


class Subframe
{
private:
  std::array<Word,10> words_;

public:
  Subframe() {}
  ~Subframe() {}
  Word& operator[](const uint8_t index);
  Word operator[](const uint8_t index) const;
  bool Bit(const uint8_t index) const;

  void Print() const;
};


class DataFrame
{
public:
  DataFrame() {}
  ~DataFrame() {}

  Subframe& operator[](uint8_t index);
  Subframe operator[](uint8_t index) const;

  void TimeIncrement(); // increases the time by one subframe
  Subframe ParityFrame(uint8_t sf); // not const because D29 and D30 are modified

  void SetSubframe(uint8_t sf_i); // zero-indexed
  void SetSubframes();
  void SetSubframe1();
  void SetSubframe2();
  void SetSubframe3();
  void SetSubframe4();
  void SetSubframe5();

  void SetTOW(uint32_t tow) { tow_ = tow; }
  void SetWeek(uint16_t week) { week_ = week; }

  // Setting other data
  void SetIntegrityFlag(const bool val) { integrity_status_flag_ = val; }
  void SetAlertFlag(const bool val) { alert_flag_ = val; }
  void SetSpoofFlag(const bool val) { anti_spoof_flag_ = val; }
  void SetL2Flag(const uint8_t val) { l2_flag_ = val; }
  void SetURA(const uint8_t val) { URA_ = val; }
  void SetHealth(const uint8_t val) { health_ = val; }
  void SetFitIntervalFlag(const bool val) { fit_interval_flag_ = val; }
  void SetAODO(const bool val) { AODO_ = val; }

  ClockData& ClockParams() { return clock_data_; }
  Ephemeris& Ephemerides() { return ephemeris_; }

  int8_t T_GD();
  uint16_t t_oc();
  int32_t a_f0();
  int16_t a_f1();
  int8_t a_f2();

  int16_t C_rs();
  int16_t del_n();
  int32_t M_0();
  int16_t C_uc();
  int32_t e();
  int16_t C_us();
  uint32_t sqrtA();
  uint16_t t_oe();
  int16_t C_ic();
  int32_t Omega_0();
  int16_t C_is();
  int32_t i_0();
  int16_t C_rc();
  int32_t omega();
  int32_t Omega_dot(); // 24 bits
  int16_t IDOT(); // 14 bits

  void RandomizeParams();

  void Print(uint8_t sf) const { subframes_[sf].Print(); }

private:
  void Preamble(const uint8_t sf_i);

  std::array<Subframe,5> subframes_;
  
  ClockData clock_data_;
  Ephemeris ephemeris_;
  uint32_t tow_; // beginning of next subframe
  uint16_t week_; // 10 bits

  bool integrity_status_flag_ = false;
  bool alert_flag_ = false;
  bool anti_spoof_flag_ = false;

  uint8_t l2_flag_ {0};
  uint8_t URA_ {0};
  uint8_t health_ {0};

  bool fit_interval_flag_ {false};
  uint8_t AODO_ {0};

  uint8_t page_ {0};
  bool D29_ {false};
  bool D30_ {false};
};

} // namespace Lnav
} // namespace Gps
#endif
