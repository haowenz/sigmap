#include "pore_model.h"

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>

#include "utils.h"

namespace sigmap {
void PoreModel::Load(const std::string &pore_model_file_path) {
  double real_start_time = GetRealTime();
  int num_kmers = 0;
  std::ifstream pore_model_file_stream(pore_model_file_path);
  std::string pore_model_file_line;
  bool first_line = true;
  while (getline(pore_model_file_stream, pore_model_file_line)) {
    std::stringstream pore_model_file_line_string_stream(pore_model_file_line);
    // skip the header
    if (pore_model_file_line[0] == '#' || pore_model_file_line.find("kmer") == 0) {
      continue;
    }
    std::string kmer;
    pore_model_file_line_string_stream >> kmer;
    if (first_line) {
      kmer_size_ = kmer.length();
      // Allocate memory to save pore model parameters
      size_t num_pore_models = 1 << (kmer_size_ * 2);
      pore_models_.assign(num_pore_models, PoreModelParameters());
      first_line = false;
    }
    assert(kmer.length() == (size_t)kmer_size_);
    uint64_t kmer_hash_value = GenerateSeedFromSequence(kmer.data(), kmer_size_, 0, kmer_size_);
    PoreModelParameters &pore_model_parameters = pore_models_[kmer_hash_value];
    pore_model_file_line_string_stream >> pore_model_parameters.level_mean >> pore_model_parameters.level_stdv >> pore_model_parameters.sd_mean >> pore_model_parameters.sd_stdv;
    ++num_kmers;
  }
  std::cerr << "Loaded " << num_kmers << " kmers in " << GetRealTime() - real_start_time << "s.\n";
}

void PoreModel::Print() {
  size_t num_pore_models = 1 << (kmer_size_ * 2);
  for (size_t i = 0; i < num_pore_models; ++i) {
    std::cerr << i << " " << pore_models_[i].level_mean << " " << pore_models_[i].level_stdv << " " << pore_models_[i].sd_mean << " " << pore_models_[i].sd_stdv << "\n";
  }
}

// This funtion return a array of level means given the [start, end) positions (0-based).
float* PoreModel::GetLevelMeansAt(const char *sequence, uint32_t start_position, uint32_t end_position) const { 
  // Note that the start and end positions should be checked before calling this function
  int32_t signal_length = end_position - start_position - kmer_size_ + 1;
  assert(signal_length > 0);
  float *signal = (float*)calloc(signal_length, sizeof(float));
  uint32_t mask = ((uint32_t)1 << (2 *kmer_size_)) - 1;
  uint32_t hash_value = GenerateSeedFromSequence(sequence, start_position + signal_length, start_position, kmer_size_);
  signal[0] = pore_models_[hash_value].level_mean;
  for (uint32_t position = start_position + 1; position < end_position - kmer_size_ + 1; ++position) {
    uint8_t current_base = CharToUint8(sequence[position + kmer_size_]); 
    if (current_base < 4) {
      hash_value = ((hash_value << 2) | current_base) & mask;
    } else {
      hash_value = (hash_value << 2) & mask;
    }
    signal[position] = pore_models_[hash_value].level_mean;
  }
  return signal;
}
} // namespace sigmap
