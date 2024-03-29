#ifndef SPATIALINDEX_H_
#define SPATIALINDEX_H_

#include <string>
#include <vector>

#include "nanoflann.hpp"
#include "sequence_batch.h"
#include "sigmap_adaptor.h"
#include "signal_batch.h"

namespace sigmap {
enum Direction {
  Negative = 0,
  Positive = 1,
};

struct SignalAnchor {
  uint32_t target_position;
  uint32_t query_position;
  float distance;
  bool operator<(const SignalAnchor &a) const {
    return std::tie(target_position, query_position, distance) <
           std::tie(a.target_position, a.query_position, a.distance);
  }
};

struct SignalAnchorChain {
  float score;
  uint32_t reference_sequence_index;
  uint32_t start_position;
  uint32_t end_position;
  uint32_t num_anchors;
  uint8_t mapq;
  Direction direction;
  // bool is_primary_chain;
  std::vector<SignalAnchor> anchors;
  bool operator>(const SignalAnchorChain &b) const {
    return std::tie(score, num_anchors, direction, reference_sequence_index,
                    start_position, end_position) >
           std::tie(b.score, b.num_anchors, b.direction,
                    b.reference_sequence_index, b.start_position,
                    b.end_position);
  }
};

class SpatialIndex {
 public:
  SpatialIndex(int min_num_seeds_required_for_mapping,
               const std::vector<int> &max_seed_frequencies,
               const std::string &index_file_path_prefix)
      : min_num_seeds_required_for_mapping_(min_num_seeds_required_for_mapping),
        max_seed_frequencies_(max_seed_frequencies),
        index_file_path_prefix_(index_file_path_prefix) {  // for read mapping
    spatial_index_ = nullptr;
  }
  SpatialIndex(int dimension, int max_leaf, int num_threads,
               const std::string &index_file_path_prefix)
      : dimension_(dimension),
        max_leaf_(max_leaf),
        num_threads_(num_threads),
        index_file_path_prefix_(
            index_file_path_prefix) {  // for index construction
    spatial_index_ = nullptr;
  }
  ~SpatialIndex() {
    if (spatial_index_ != nullptr) {
      delete spatial_index_;
    }
  }
  void GeneratePointCloudOnOneDirection(
      Direction direction, uint32_t signal_index,
      const std::vector<std::vector<bool> > &is_masked,
      const float *signal_values, size_t signal_length, int step_size,
      std::vector<Point> &point_cloud);
  // void GetSignalIndexAndPosition(size_t point_index, size_t num_signals,
  // const std::vector<std::vector<float> > &signals, size_t &signal_index,
  // size_t &signal_position);
  void Construct(size_t num_signals,
                 const std::vector<std::vector<bool> > &positive_is_masked,
                 const std::vector<std::vector<bool> > &negative_is_masked,
                 const std::vector<std::vector<float> > &positive_signals,
                 const std::vector<std::vector<float> > &negative_signals);
  void Save();
  void Load();

  void GenerateChains(const std::vector<float> &query_signal,
                      const std::vector<float> &query_signal_stdvs,
                      uint32_t query_start_offset,
                      int query_point_cloud_step_size, float search_radius,
                      size_t num_target_signals,
                      std::vector<SignalAnchorChain> &chains);
  void TracebackChains(
      int min_num_anchors, Direction direction, size_t chain_end_anchor_index,
      uint32_t chain_target_signal_index,
      const std::vector<float> &chaining_scores,
      const std::vector<size_t> &chaining_predecessors,
      const std::vector<std::vector<SignalAnchor> > &anchors_on_diff_signals,
      std::vector<bool> &anchor_is_used,
      std::vector<SignalAnchorChain> &chains);
  void GeneratePrimaryChains(std::vector<SignalAnchorChain> &chains);
  void ComputeMAPQ(std::vector<SignalAnchorChain> &chains);

 protected:
  int dimension_;
  int max_leaf_ = 10;
  // int window_size_;
  int min_num_seeds_required_for_mapping_;
  std::vector<int> max_seed_frequencies_;
  int num_threads_;
  std::string index_file_path_prefix_;  // For now, I have to use two files *.pt
                                        // for dim, max leaf, points and *.si
                                        // for spatial index.
  SigmapAdaptor<float> *spatial_index_;
  std::vector<Point> point_cloud_;  // TODO(Haowen): remove this duplicate later
};
}  // namespace sigmap

#endif  // SPATIALINDEX_H_
