#include "spatial_index.h"

#include <assert.h>

#include <algorithm>
#include <iostream>

#include "utils.h"

namespace sigmap {
bool compare(const std::pair<float, size_t> &left,
             const std::pair<float, size_t> &right) {
  if (left.first > right.first) {
    return true;
  } else if (left.first == right.first) {
    return (left.second > right.second);
  } else {
    return false;
  }
}

bool compare1(const std::pair<float, size_t> &left,
              const std::pair<float, size_t> &right) {
  if (left.first > right.first) {
    return true;
  } else if (left.first == right.first) {
    return (left.second < right.second);
  } else {
    return false;
  }
}

void SpatialIndex::GeneratePointCloudOnOneDirection(
    Direction direction, uint32_t signal_index,
    const std::vector<std::vector<bool> > &is_masked,
    const float *signal_values, size_t signal_length, int step_size,
    std::vector<Point> &point_cloud) {
  if (signal_length >= (uint32_t)dimension_) {
    for (size_t signal_position = 0;
         signal_position < signal_length - dimension_ + 1;
         signal_position += step_size) {
      if (!is_masked[signal_index][signal_position]) {
        if (signal_position == 0 || point_cloud.empty() ||
            (signal_position > 0 && abs(signal_values[signal_position] -
                                        point_cloud.back().value) > 0.01)) {
          uint64_t strand = direction == Positive ? 0 : 1;
          uint64_t position =
              ((((uint64_t)signal_index) << 32 | (uint32_t)signal_position)
               << 1) |
              strand;
          point_cloud.emplace_back(position, signal_values[signal_position]);
        }
      }
    }
  }
}

// void SpatialIndex::GetSignalIndexAndPosition(size_t point_index, size_t
// num_signals, const std::vector<std::vector<float> > &signals, size_t
// &signal_index, size_t &signal_position) {
//  for (signal_index = 0; signal_index < num_signals; ++signal_index) {
//    signal_position = signals[signal_index].size() - dimension_ + 1;
//    if (point_index > signal_position) {
//      point_index -= signal_position;
//    } else {
//      signal_position = point_index;
//      break;
//    }
//  }
//}

void SpatialIndex::Construct(
    size_t num_signals,
    const std::vector<std::vector<bool> > &positive_is_masked,
    const std::vector<std::vector<bool> > &negative_is_masked,
    const std::vector<std::vector<float> > &positive_signals,
    const std::vector<std::vector<float> > &negative_signals) {
  double real_start_time = GetRealTime();
  int signal_point_step_size = 1;
  point_cloud_.reserve(2 * GetSignalsTotalLength(positive_signals));
  for (size_t signal_index = 0; signal_index < num_signals; ++signal_index) {
    GeneratePointCloudOnOneDirection(Positive, signal_index, positive_is_masked,
                                     positive_signals[signal_index].data(),
                                     positive_signals[signal_index].size(),
                                     signal_point_step_size, point_cloud_);
  }
  for (size_t signal_index = 0; signal_index < num_signals; ++signal_index) {
    GeneratePointCloudOnOneDirection(Negative, signal_index, negative_is_masked,
                                     negative_signals[signal_index].data(),
                                     negative_signals[signal_index].size(),
                                     signal_point_step_size, point_cloud_);
  }
  std::cerr << "Collected " << point_cloud_.size() << " points.\n";
  // sort(point_cloud_.begin(), point_cloud_.end());
  // std::cerr << "Sorted " << point_cloud_.size() << " points.\n";
  spatial_index_ = new SigmapAdaptor<float>(dimension_ /*dim*/, point_cloud_,
                                            max_leaf_ /* max leaf */);
  spatial_index_->index->buildIndex();
  std::cerr << "Built spatial index.\n";
  std::cerr << "Built index successfully in " << GetRealTime() - real_start_time
            << "s.\n";
}

void SpatialIndex::Save() {
  double real_start_time = GetRealTime();
  FILE *point_cloud_file =
      fopen((index_file_path_prefix_ + ".pt").c_str(), "wb");
  assert(point_cloud_file != NULL);
  fwrite(&dimension_, sizeof(int), 1, point_cloud_file);
  fwrite(&max_leaf_, sizeof(int), 1, point_cloud_file);
  size_t point_cloud_size = point_cloud_.size();
  fwrite(&point_cloud_size, sizeof(size_t), 1, point_cloud_file);
  // TODO(Haowen): we don't have to save the whole point cloud. Instead, we can
  // just save the feature signal and build the point cloud on the fly. I will
  // leave this as an easy later code refactor work.
  // for (size_t pi = 0; pi < point_cloud_size; ++pi) {
  //  fwrite(point_cloud_[pi].data(), sizeof(float), dimension_,
  //  point_cloud_file);
  //}
  fwrite(point_cloud_.data(), sizeof(Point), point_cloud_size,
         point_cloud_file);
  fclose(point_cloud_file);
  FILE *spatial_index_file =
      fopen((index_file_path_prefix_ + ".si").c_str(), "wb");
  assert(spatial_index_file != NULL);
  spatial_index_->index->saveIndex(spatial_index_file);
  fclose(spatial_index_file);
  std::cerr << "Saved in " << GetRealTime() - real_start_time << "s.\n";
}

void SpatialIndex::Load() {
  double real_start_time = GetRealTime();
  FILE *point_cloud_file =
      fopen((index_file_path_prefix_ + ".pt").c_str(), "rb");
  fread(&dimension_, sizeof(int), 1, point_cloud_file);
  fread(&max_leaf_, sizeof(int), 1, point_cloud_file);
  size_t point_cloud_size = 0;
  fread(&point_cloud_size, sizeof(size_t), 1, point_cloud_file);
  point_cloud_.resize(point_cloud_size);
  // for (size_t pi = 0; pi < point_cloud_size; ++pi) {
  //  point_cloud_[pi].resize(dimension_);
  //  fread(point_cloud_[pi].data(), sizeof(float), dimension_,
  //  point_cloud_file);
  //}
  fread(point_cloud_.data(), sizeof(Point), point_cloud_size, point_cloud_file);
  fclose(point_cloud_file);
  std::cerr << "Load point cloud successfully! dim: " << dimension_
            << ", max leaf: " << max_leaf_
            << ", point cloud size: " << point_cloud_size << "\n";
  FILE *spatial_index_file =
      fopen((index_file_path_prefix_ + ".si").c_str(), "rb");
  assert(spatial_index_file != NULL);
  spatial_index_ = new SigmapAdaptor<float>(dimension_ /*dim*/, point_cloud_,
                                            max_leaf_ /* max leaf */);
  spatial_index_->index->loadIndex(spatial_index_file);
  size_t used_memory =
      spatial_index_->index->usedMemory(*(spatial_index_->index));
  std::cerr << "Memory: " << used_memory << ".\n";
  fclose(spatial_index_file);
  std::cerr << "Loaded index successfully in "
            << GetRealTime() - real_start_time << "s.\n";
}

void SpatialIndex::TracebackChains(
    int min_num_anchors, Direction direction, size_t chain_end_anchor_index,
    uint32_t chain_target_signal_index,
    const std::vector<float> &chaining_scores,
    const std::vector<size_t> &chaining_predecessors,
    const std::vector<std::vector<SignalAnchor> > &anchors_on_diff_signals,
    std::vector<bool> &anchor_is_used, std::vector<SignalAnchorChain> &chains) {
  if (!anchor_is_used[chain_end_anchor_index]) {
    std::vector<SignalAnchor> anchors;
    anchors.reserve(100);
    bool stop_at_an_used_anchor = false;
    size_t chain_start_anchor_index = chain_end_anchor_index;
    // Add the end anchor
    anchors.push_back(anchors_on_diff_signals[chain_target_signal_index]
                                             [chain_start_anchor_index]);
    // The end anchor has an used predecessor
    if (chaining_predecessors[chain_start_anchor_index] !=
            chain_start_anchor_index &&
        anchor_is_used[chaining_predecessors[chain_start_anchor_index]]) {
      stop_at_an_used_anchor = true;
    }
    anchor_is_used[chain_start_anchor_index] = true;
    uint32_t chain_num_anchors = 1;
    while (chaining_predecessors[chain_start_anchor_index] !=
               chain_start_anchor_index &&
           !anchor_is_used[chaining_predecessors[chain_start_anchor_index]]) {
      chain_start_anchor_index =
          chaining_predecessors[chain_start_anchor_index];
      anchors.push_back(anchors_on_diff_signals[chain_target_signal_index]
                                               [chain_start_anchor_index]);
      if (chaining_predecessors[chain_start_anchor_index] !=
              chain_start_anchor_index &&
          anchor_is_used[chaining_predecessors[chain_start_anchor_index]]) {
        stop_at_an_used_anchor = true;
      }
      anchor_is_used[chain_start_anchor_index] = true;
      ++chain_num_anchors;
    }
    if (chain_num_anchors >= (uint32_t)min_num_anchors) {
      float adjusted_chaining_score = chaining_scores[chain_end_anchor_index];
      if (stop_at_an_used_anchor) {
        adjusted_chaining_score -=
            chaining_scores[chaining_predecessors[chain_start_anchor_index]];
      }
      chains.emplace_back(
          SignalAnchorChain{adjusted_chaining_score, chain_target_signal_index,
                            anchors_on_diff_signals[chain_target_signal_index]
                                                   [chain_start_anchor_index]
                                                       .target_position,
                            anchors_on_diff_signals[chain_target_signal_index]
                                                   [chain_end_anchor_index]
                                                       .target_position,
                            chain_num_anchors, 0, direction, anchors});
    }
  }
}

void SpatialIndex::GeneratePrimaryChains(
    std::vector<SignalAnchorChain> &chains) {
  std::sort(chains.begin(), chains.end(), std::greater<SignalAnchorChain>());
  std::vector<SignalAnchorChain> primary_chains;
  primary_chains.reserve(chains.size());
  primary_chains.emplace_back(chains[0]);
  for (uint32_t ci = 1; ci < chains.size(); ++ci) {
    bool is_primary = true;
    if (chains[ci].score < primary_chains.back().score / 3) {
      break;
    } else {
      for (uint32_t pi = 0; pi < primary_chains.size(); ++pi) {
        if (chains[ci].reference_sequence_index ==
            primary_chains[pi].reference_sequence_index) {
          if (std::max(chains[ci].start_position,
                       primary_chains[pi].start_position) >
              std::min(chains[ci].end_position,
                       primary_chains[pi].end_position)) {
            // primary_chains[pi].score += chains[ci].score;
          } else {
            is_primary = false;
            break;
          }
        }
      }
    }
    if (is_primary) {
      primary_chains.emplace_back(chains[ci]);
    }
  }
  chains.swap(primary_chains);
}

void SpatialIndex::ComputeMAPQ(std::vector<SignalAnchorChain> &chains) {
  if (chains.size() == 1) {
    chains[0].mapq = 60;
    return;
  } else {
    int mapq =
        40 *
        (1 -
         chains[1].score /
             chains[0].score);  // * std::min((size_t)1, chains[0].num_anchors /
                                // 20) * log(chains[0].score);
    if (mapq > 60) {
      mapq = 60;
    }
    if (mapq < 0) {
      mapq = 0;
    }
    chains[0].mapq = (uint8_t)mapq;
  }
}

void SpatialIndex::GenerateChains(const std::vector<float> &query_signal,
                                  const std::vector<float> &query_signal_stdvs,
                                  uint32_t query_start_offset,
                                  int query_point_cloud_step_size,
                                  float search_radius,
                                  size_t num_target_signals,
                                  std::vector<SignalAnchorChain> &chains) {
  // Chaining parameters
  int max_gap_length = 2000;
  int max_target_gap_length = 5000;
  int chaining_band_length = 5000;
  int max_num_skips = 25;
  int min_num_anchors = 2;
  int num_best_chains = 3;
  int num_nearest_points = 5000;
  float min_chaining_score = 10;
  std::vector<std::vector<std::vector<SignalAnchor> > > anchors_on_diff_signals(
      2);
  anchors_on_diff_signals[0] =
      std::vector<std::vector<SignalAnchor> >(num_target_signals);
  anchors_on_diff_signals[1] =
      std::vector<std::vector<SignalAnchor> >(num_target_signals);
  std::vector<std::vector<SignalAnchor> > &positive_anchors_on_diff_signals =
      anchors_on_diff_signals[0];
  std::vector<std::vector<SignalAnchor> > &negative_anchors_on_diff_signals =
      anchors_on_diff_signals[1];
  // Get anchors in previous chains
  std::vector<SignalAnchorChain> previous_chains;
  previous_chains.swap(chains);
  if (previous_chains.size() > 0) {
    for (uint32_t chain_index = 0; chain_index < previous_chains.size();
         ++chain_index) {
      int strand = previous_chains[chain_index].direction == Positive ? 0 : 1;
      uint32_t reference_sequence_index =
          previous_chains[chain_index].reference_sequence_index;
      assert(previous_chains[chain_index].num_anchors ==
             previous_chains[chain_index].anchors.size());
      anchors_on_diff_signals[strand][reference_sequence_index].reserve(
          previous_chains[chain_index].num_anchors);
      for (uint32_t anchor_index = 0;
           anchor_index < previous_chains[chain_index].num_anchors;
           ++anchor_index) {
        anchors_on_diff_signals[strand][reference_sequence_index].emplace_back(
            previous_chains[chain_index].anchors[anchor_index]);
      }
    }
  }
  nanoflann::SearchParams params;
  params.sorted = false;
  std::vector<std::pair<size_t, float> > point_anchors;
  // Find reliable seeds
  std::vector<std::pair<float, size_t> > mean_diff_position;
  mean_diff_position.reserve(
      (query_signal.size() - dimension_ + 1) / query_point_cloud_step_size + 1);
  for (uint32_t pi = 0; pi < query_signal.size() - dimension_ + 1; ++pi) {
    float min_diff = std::numeric_limits<float>::max();
    // for (int di = 0; di < dimension_; ++di) {
    for (int di = 1; di < dimension_; ++di) {
      float diff = abs(query_signal[pi + di] - query_signal[pi + di - 1]);
      // if (diff < min_diff) {
      min_diff += diff;
      //}
      // float diff = query_signal_stdvs[pi + di];
      // if (diff < min_diff) {
      //  min_diff = diff;
      //}
      // min_diff += query_signal_stdvs[pi + di];
    }
    mean_diff_position.emplace_back(min_diff, pi);
  }
  std::sort(mean_diff_position.begin(), mean_diff_position.end(), compare1);
  // std::sort(mean_diff_position.begin(), mean_diff_position.end());
  // Collect anchors
  uint32_t previous_position = 0;
  uint32_t num_positions = 0;
  for (uint32_t pi = 0; pi < query_signal.size() - dimension_ + 1; ++pi) {
    uint32_t position = mean_diff_position[pi].second;
    if (position < previous_position + query_point_cloud_step_size &&
        position + query_point_cloud_step_size > previous_position) {
      continue;
    }
    // std::cout << position << "\n";
    // size_t num_results = 100;
    // std::vector<size_t> ret_indexes(num_results);
    // std::vector<float> out_dists_sqr(num_results);
    // nanoflann::KNNResultSet<float> resultSet(num_results);
    // resultSet.init(&ret_indexes[0], &out_dists_sqr[0] );
    // size_t num_point_anchors =
    // spatial_index_->index->findNeighbors(resultSet, query_signal.data() +
    // position, nanoflann::SearchParams(10));
    size_t num_point_anchors = spatial_index_->index->radiusSearch(
        query_signal.data() + position, search_radius, point_anchors, params);
    // std::cout << "radiusSearch(): radius=" << search_radius << " -> " <<
    // num_point_anchors << " matches\n"; for (size_t ai = 0; ai < num_results;
    // ai++) { float previous_distance = dimension_;
    for (size_t ai = 0;
         ai < num_point_anchors && ai < (uint32_t)num_nearest_points; ai++) {
      // if (point_anchors[ai].second > previous_distance * 2) {
      //  break;
      //}
      // std::cout << "cloud size =" << point_cloud_.size() << " -> " <<
      // point_anchors[ai].first << " matches\n";
      Point &point = point_cloud_[point_anchors[ai].first];
      // Point &point = point_cloud_[ret_indexes[ai]];
      uint32_t target_signal_index = point.position >> 33,
               target_signal_position = point.position >> 1;
      Direction target_signal_direction =
          (point.position & 1) == 0 ? Positive : Negative;
      // GetSignalIndexAndPosition(point_anchors[ai].first, num_target_signals,
      // target_signals, target_signal_index, target_signal_position);
      if (target_signal_direction == Positive) {
        positive_anchors_on_diff_signals[target_signal_index].emplace_back(
            SignalAnchor{target_signal_position, position + query_start_offset,
                         point_anchors[ai].second});
        // positive_anchors_on_diff_signals[target_signal_index].emplace_back(SignalAnchor{target_signal_position,
        // position, out_dists_sqr[ai]});
      } else {
        negative_anchors_on_diff_signals[target_signal_index].emplace_back(
            SignalAnchor{target_signal_position, position + query_start_offset,
                         point_anchors[ai].second});
        // negative_anchors_on_diff_signals[target_signal_index].emplace_back(SignalAnchor{target_signal_position,
        // position, out_dists_sqr[ai]});
      }
      // previous_distance = point_anchors[ai].second;
      // std::cout << "idx["<< ai << "]=" << point_anchors[ai].first << "
      // dist["<< ai << "]=" << point_anchors[ai].second << std::endl;
    }
    ++num_positions;
    if (num_positions >=
        (query_signal.size() - dimension_ + 1) / query_point_cloud_step_size) {
      break;
    }
    previous_position = position;
  }
  // Sort the anchors based on their occurrence on target signal
  for (size_t target_signal_index = 0; target_signal_index < num_target_signals;
       ++target_signal_index) {
    std::sort(positive_anchors_on_diff_signals[target_signal_index].begin(),
              positive_anchors_on_diff_signals[target_signal_index].end());
    std::sort(negative_anchors_on_diff_signals[target_signal_index].begin(),
              negative_anchors_on_diff_signals[target_signal_index].end());
  }
  // Chaining DP done on each individual target signal
  float max_chaining_score = 0;  // std::numeric_limits<float>::min();
  for (size_t target_signal_index = 0; target_signal_index < num_target_signals;
       ++target_signal_index) {
    for (int direction_i = 0; direction_i < 2; ++direction_i) {
      std::vector<float> chaining_scores;
      chaining_scores.reserve(
          anchors_on_diff_signals[direction_i][target_signal_index].size());
      std::vector<size_t> chaining_predecessors;
      chaining_predecessors.reserve(
          anchors_on_diff_signals[direction_i][target_signal_index].size());
      std::vector<bool> anchor_is_used;
      anchor_is_used.reserve(
          anchors_on_diff_signals[direction_i][target_signal_index].size());
      std::vector<std::pair<float, size_t> > end_anchor_index_chaining_scores;
      end_anchor_index_chaining_scores.reserve(10);
      for (size_t anchor_index = 0;
           anchor_index <
           anchors_on_diff_signals[direction_i][target_signal_index].size();
           ++anchor_index) {
        float distance_coefficient =
            1 - 0.2 *
                    anchors_on_diff_signals[direction_i][target_signal_index]
                                           [anchor_index]
                                               .distance /
                    search_radius;
        chaining_scores.emplace_back(distance_coefficient * dimension_);
        chaining_predecessors.emplace_back(anchor_index);
        anchor_is_used.push_back(false);
        int32_t current_anchor_target_position =
            anchors_on_diff_signals[direction_i][target_signal_index]
                                   [anchor_index]
                                       .target_position;
        int32_t current_anchor_query_position =
            anchors_on_diff_signals[direction_i][target_signal_index]
                                   [anchor_index]
                                       .query_position;
        int32_t start_anchor_index = 0;
        if (anchor_index > (size_t)chaining_band_length) {
          start_anchor_index = anchor_index - chaining_band_length;
        }
        int32_t previous_anchor_index = anchor_index - 1;
        int32_t num_skips = 0;
        for (; previous_anchor_index >= start_anchor_index;
             --previous_anchor_index) {
          // std::cerr << previous_anchor_index << " " << start_anchor_index <<
          // " " << anchor_index << "\n";
          int32_t previous_anchor_target_position =
              anchors_on_diff_signals[direction_i][target_signal_index]
                                     [previous_anchor_index]
                                         .target_position;
          int32_t previous_anchor_query_position =
              anchors_on_diff_signals[direction_i][target_signal_index]
                                     [previous_anchor_index]
                                         .query_position;
          if (previous_anchor_query_position == current_anchor_query_position) {
            continue;
          }
          if (previous_anchor_target_position ==
              current_anchor_target_position) {
            continue;
          }
          if (previous_anchor_target_position + max_target_gap_length <
              current_anchor_target_position) {
            break;
          }
          int32_t target_position_diff =
              current_anchor_target_position - previous_anchor_target_position;
          assert(target_position_diff > 0);
          int32_t query_position_diff =
              current_anchor_query_position - previous_anchor_query_position;
          float current_chaining_score = 0;
          // std::cerr << "curr_ai: " << anchor_index << ", adjusted_d: " <<
          // dimension_ * distance_coefficient << ", pre_ai: " <<
          // previous_anchor_index << ", "; std::cerr << "target_diff: " <<
          // target_position_diff << ", query_diff: " << query_position_diff <<
          // "\n";
          if (query_position_diff < 0) {
            // current_chaining_score = std::numeric_limits<float>::min();
            continue;
          } else {
            float matching_dimensions =
                std::min(std::min(target_position_diff, query_position_diff),
                         dimension_) *
                distance_coefficient;
            int gap_length =
                std::abs(target_position_diff - query_position_diff);
            float gap_scale =
                target_position_diff > 0
                    ? (float)query_position_diff / target_position_diff
                    : 1;
            // float gap_cost = 0;
            // if (gap_length != 0) {
            if (gap_length < max_gap_length && gap_scale < 5 &&
                gap_scale > 0.75) {
              current_chaining_score = chaining_scores[previous_anchor_index] +
                                       matching_dimensions;  // - gap_cost;
            } else {
              // gap_cost = std::numeric_limits<float>::max();
              // current_chaining_score =
              // std::numeric_limits<float>::min();//chaining_scores[previous_anchor_index]
              // + matching_dimensions - gap_cost;
            }
            //}
            // std::cerr << ", matching_d: " << matching_dimensions;
            // std::cerr << ", gap_len: " << gap_length;
            // std::cerr << ", gap_cost: " << gap_cost;
            // std::cerr << ", current chaining score: " <<
            // current_chaining_score << "\n";
          }
          if (current_chaining_score > chaining_scores[anchor_index]) {
            chaining_scores[anchor_index] = current_chaining_score;
            chaining_predecessors[anchor_index] = previous_anchor_index;
            --num_skips;
            // std::cerr << "update_curr_max: " << current_chaining_score <<
            // "\n";
          } else {
            ++num_skips;
            if (num_skips > max_num_skips) {
              break;
            }
          }
        }
        // Update chain with max score
        if (chaining_scores[anchor_index] > max_chaining_score) {
          max_chaining_score = chaining_scores[anchor_index];
        }
        if (chaining_scores.back() >= min_chaining_score &&
            chaining_scores.back() > max_chaining_score / 2) {
          end_anchor_index_chaining_scores.emplace_back(chaining_scores.back(),
                                                        anchor_index);
        }
      }
      // Sort a vector of <end anchor index, chaining score>
      std::sort(end_anchor_index_chaining_scores.begin(),
                end_anchor_index_chaining_scores.end(), compare);
      // Traceback all chains from higest score to lowest
      for (size_t anchor_index = 0;
           anchor_index < end_anchor_index_chaining_scores.size() &&
           anchor_index < (size_t)num_best_chains;
           ++anchor_index) {
        TracebackChains(
            min_num_anchors, direction_i == 0 ? Positive : Negative,
            end_anchor_index_chaining_scores[anchor_index].second,
            target_signal_index, chaining_scores, chaining_predecessors,
            anchors_on_diff_signals[direction_i], anchor_is_used, chains);
        if (chaining_scores[end_anchor_index_chaining_scores[anchor_index]
                                .second] < max_chaining_score / 2) {
          break;
        }
      }
    }
  }
  if (chains.size() > 0) {
    // Generate primary chains
    GeneratePrimaryChains(chains);
    // Compute MAPQ
    ComputeMAPQ(chains);
  }
}
}  // namespace sigmap
