#ifndef SIGMAP_H_
#define SIGMAP_H_

#include <memory>
#include <string>

#include "event.h"
#include "output_tools.h"
#include "signal_batch.h"
#include "utils.h"

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
  void ConstructIndex();
  // for mapping
  Sigmap(int num_threads, const std::string &reference_file_path, const std::string &pore_model_file_path, const std::string &signal_directory, const std::string &reference_index_file_path, const std::string &output_file_path) : num_threads_(num_threads), reference_file_path_(reference_file_path), pore_model_file_path_(pore_model_file_path), signal_directory_(signal_directory), reference_index_file_path_(reference_index_file_path), output_file_path_(output_file_path) {}
  Sigmap(float search_radius, int read_seeding_step_size, int num_threads, int max_num_chunks, int stop_mapping_min_num_anchors, int output_mapping_min_num_anchors, float stop_mapping_ratio, float output_mapping_ratio, float stop_mapping_mean_ratio, float output_mapping_mean_ratio, const std::string &reference_file_path, const std::string &pore_model_file_path, const std::string &signal_directory, const std::string &reference_index_file_path, const std::string &output_file_path) : search_radius_(search_radius), read_seeding_step_size_(read_seeding_step_size), num_threads_(num_threads), max_num_chunks_(max_num_chunks), stop_mapping_min_num_anchors_(stop_mapping_min_num_anchors), output_mapping_min_num_anchors_(output_mapping_min_num_anchors), stop_mapping_ratio_(stop_mapping_ratio), output_mapping_ratio_(output_mapping_ratio), stop_mapping_mean_ratio_(stop_mapping_mean_ratio), output_mapping_mean_ratio_(output_mapping_mean_ratio), reference_file_path_(reference_file_path), pore_model_file_path_(pore_model_file_path), signal_directory_(signal_directory), reference_index_file_path_(reference_index_file_path), output_file_path_(output_file_path) {}
  void Map();
  void StreamingMap();
  void DTWAlign();
  void CWTAlign();
  // Support functions
  void GenerateEvents(size_t start, size_t end, const Signal &signal, std::vector<float> &feature_signal, std::vector<float> &feature_signal_stdvs);
  void GenerateFeatureSignalUsingCWT(const Signal &signal, float scale0, std::vector<float> &feature_signal, std::vector<size_t> &feature_positions);
  float GenerateMADNormalizedSignal(const float *signal_values, size_t signal_length, std::vector<float> &normalized_signal);
  float GenerateZscoreNormalizedSignal(const float *signal_values, size_t signal_length, std::vector<float> &normalized_signal);
  void GenerateCWTSignal(const float *signal_values, size_t signal_length, float scale0, std::vector<float> &cwt_signal);
  void GeneratePeaks(const float *signal_values, size_t signal_length, float selective, std::vector<float> &peaks, std::vector<size_t> &peak_positions);
  float sDTW(const Signal &target_signal, const Signal &query_signal);
  float sDTW(const float *target_signal_values, size_t target_length, const float *query_signal_values, size_t query_length, ssize_t &mapping_end_position);
  //void EmplaceBackMappingRecord(uint32_t read_id, const char *read_name, uint32_t read_length, uint32_t barcode, uint32_t fragment_start_position, uint32_t fragment_length, uint8_t mapq, uint8_t direction, uint8_t is_unique, std::vector<PAFMapping> *mappings_on_diff_ref_seqs);
  void OutputMappingsInVector(uint8_t mapq_threshold, uint32_t num_reference_sequences, const SequenceBatch &reference, const std::vector<std::vector<PAFMapping> > &mappings);
  uint32_t MoveMappingsInBuffersToMappingContainer(uint32_t num_reference_sequences, std::vector<std::vector<std::vector<PAFMapping> > > *mappings_on_diff_ref_seqs_for_diff_threads_for_saving);
  void GenerateMaskedPositions(int kmer_size, float frequency, uint32_t num_reference_sequences, const SequenceBatch &sequence_batch, std::vector<std::vector<bool> > &positive_is_masked, std::vector<std::vector<bool> > &negative_is_masked);
  // Output debug files
  void FAST5ToText();
  void EventsToText();

 protected:
  int dimension_;
  int max_leaf_;
  float search_radius_;
  int read_seeding_step_size_;
  int num_threads_;
  int max_num_chunks_;
  int stop_mapping_min_num_anchors_;
  int output_mapping_min_num_anchors_;
  float stop_mapping_ratio_;
  float output_mapping_ratio_;
  float stop_mapping_mean_ratio_;
  float output_mapping_mean_ratio_;
  std::string reference_file_path_;
  std::string pore_model_file_path_;
  std::string signal_directory_;
  std::string reference_index_file_path_;
  std::string output_file_path_;
  std::unique_ptr<OutputTools<PAFMapping> > output_tools_;
  std::vector<std::vector<PAFMapping> > mappings_on_diff_ref_seqs_;
};
} // namespace sigmap

#endif // SIGMAP_H_
