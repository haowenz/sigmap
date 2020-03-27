#ifndef SIGMAP_H_
#define SIGMAP_H_

#include <string>

#include "signal_batch.h"

namespace sigmap {
class SigmapDriver {
 public:
  SigmapDriver() {}
  ~SigmapDriver() {}
  void ParseArgsAndRun(int argc, char *argv[]);
};

class Sigmap {
 public:
  Sigmap(){}
  ~Sigmap(){}
  // for index construction
  Sigmap(int dimension, int max_leaf, const std::string &reference_file_path, const std::string &pore_model_file_path, const std::string &output_file_path) : dimension_(dimension), max_leaf_(max_leaf), reference_file_path_(reference_file_path), pore_model_file_path_(pore_model_file_path), output_file_path_(output_file_path) {}
  // for mapping
  Sigmap(const std::string &reference_file_path, const std::string &pore_model_file_path, const std::string &signal_directory, const std::string &reference_index_file_path, const std::string &output_file_path) : reference_file_path_(reference_file_path), pore_model_file_path_(pore_model_file_path), signal_directory_(signal_directory), reference_index_file_path_(reference_index_file_path), output_file_path_(output_file_path) {}
  void Map();
  void DTWAlign();
  void CWTAlign();
  void ConstructIndex();
  void GetFeatureSignal(const Signal &signal, float scale0, std::vector<float> &feature_signal, std::vector<size_t> &feature_positions);
  void FAST5ToText();
  float sDTW(const Signal &target_signal, const Signal &query_signal);
  float sDTW(const float *target_signal, size_t target_length, const float *query_signal, size_t query_length, ssize_t &mapping_end_position);
 protected:
  int dimension_;
  int max_leaf_;
  std::string reference_file_path_;
  std::string pore_model_file_path_;
  std::string signal_directory_;
  std::string reference_index_file_path_;
  std::string output_file_path_;
};
} // namespace sigmap

#endif // SIGMAP_H_
