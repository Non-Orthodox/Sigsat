#include <iostream>
#include <array>
#include <vector>
#include <fstream>
#include <complex>

#include <Eigen/Dense>

#include "gps_common.hpp"
#include "python_plotting.hpp"

/*
The purpose of this test is to ensure that numerical errors do not significantly affect Baseband CA code sampling.
Additionally, this tests to ensure that no code phase drift is present across multiple sampled intervals.
This test demonstrates that numerical errors in code phase do not compound over time.
Although integer-millisecond intervals are used for simple verification, the principle should apply to any arbitrary interval.
*/
void BasebandCaSamplingTest()
{
  double chip_tolerance = 1.0e-10; // corresponds to approx 1.0e-16 seconds
  std::array<bool,1023> ca_code;
  Gps::GenCA(&ca_code, 1);

  std::cout << "Baseband CA Sampling Test: ";
  bool passed = true;
  std::size_t num_chip_starts = 97; // prime number so that precision errors must exist for every test
  for (std::size_t i = 0; i < num_chip_starts; i++) { 
    double chip = static_cast<double>(i) / static_cast<double>(num_chip_starts);
    double init_chip = chip;
    int8_t samples [1000000];
    for (std::size_t j = 0; j < 10; j++) {
      Gps::SampleBasebandCa<int8_t>(ca_code, samples, 1000000, 10.0e6, chip, 100);
      std::cout << (chip - init_chip) << '\n';
    }
    double diff = chip - init_chip;
    passed &= ( (diff > 0.0 ? diff : -diff) < chip_tolerance );
    if (!passed) break;
  }
  if (passed) std::cout << "passed\n";
  else std::cout << "failed\n";
}


/*
This test verifies the expected properties of correlating CA codes with themselves
for various starting chip values. This only tests baseband CA code data with no doppler. 
*/
void CaCorrelationTest(const double f_s, double true_chip, const uint8_t prn)
{
  true_chip = fmod(true_chip,1023.0);
  std::size_t arr_size = static_cast<std::size_t>((f_s / 1000.0) + 0.5);
  std::array<bool,1023> ca_code;
  Gps::GenCA(&ca_code, prn);
  std::vector<double> samples(arr_size);
  Gps::SampleBasebandCa<double>(ca_code, samples.data(), arr_size, f_s, true_chip, 1.0);
  
  std::size_t num_chip_starts = 1023 * 2;
  std::vector<double> chips(num_chip_starts);
  std::vector<double> corr_outputs(num_chip_starts);
  for (std::size_t i = 0; i < num_chip_starts; i++) {
    chips[i] = 1023.0 * static_cast<double>(i) / static_cast<double>(num_chip_starts);
    std::vector<double> replica(arr_size);
    Gps::SampleBasebandCa<double>(ca_code, replica.data(), replica.size(), f_s, chips[i], 1.0);
    corr_outputs[i] = Gps::Correlate(samples, replica);
  }

  PythonPlot plt;
  plt.plot(chips,corr_outputs,"test plot");
  plt.show();
  plt.execute();
}


void ComplexCaCorrelationTest()
{
  double f_s = 10.0e6;
  double true_chip = 100.0;
  double true_carrier_phase = 30.0 * Gps::PI / 180.0;
  uint8_t prn = 1;
  
  true_chip = fmod(true_chip,1023.0);
  std::size_t arr_size = static_cast<std::size_t>((f_s / 1000.0) + 0.5);
  std::array<bool,1023> ca_code;
  Gps::GenCA(&ca_code, prn);
  std::vector<std::complex<double>> samples(arr_size);
  Gps::SampleCa(ca_code, samples.data(), arr_size, f_s, Gps::CA_RATE, true_chip,
                0.0, true_carrier_phase, 1.0);
  
  std::size_t num_chip_starts = 1023 * 10;
  std::vector<double> chips(num_chip_starts);
  std::vector<double> corr_outputs(num_chip_starts);
  for (std::size_t i = 0; i < num_chip_starts; i++) {
    chips[i] = 1023.0 * static_cast<double>(i) / static_cast<double>(num_chip_starts);
    std::vector<double> replica(arr_size);
    Gps::SampleBasebandCa(ca_code, replica.data(), replica.size(), f_s, chips[i], 1.0);
    corr_outputs[i] = sqrt(norm(Gps::ComplexCorrelate(samples, replica)));
  }

  PythonPlot plt;
  plt.plot(chips,corr_outputs,"test plot");
  plt.show();
  plt.execute();
}


int main()
{
  BasebandCaSamplingTest();
  CaCorrelationTest(10.0e6, 100.0, 1);
  ComplexCaCorrelationTest();
  return 0;
}
/*
Next to add:
  Sample baseband data + code
  general FFT stuff


  Data Frame class:
    -contains all relevent data for full frame generation (almanac optional)
    -can generate subframes that can be interacted with like an array of boolean values
    -can be updated with individual parameters OR raw subframe data
    -
*/

