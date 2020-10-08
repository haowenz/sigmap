#include "spatial_index.h"

#include <algorithm>
#include <assert.h>
#include <iostream>

#include "utils.h"

namespace sigmap {
bool compare(const std::pair<float,size_t> &left, const std::pair<float,size_t> &right) {
  if (left.first > right.first) {
    return true;
  } else if (left.first == right.first) {
    return (left.second > right.second);
  } else { 
    return false;
  }
}

bool compare1(const std::pair<float,size_t> &left, const std::pair<float,size_t> &right) {
  if (left.first > right.first) {
    return true;
  } else if (left.first == right.first) {
    return (left.second < right.second);
  } else { 
    return false;
  }
}
//void Index::Statistics(uint32_t num_sequences, const SequenceBatch &reference) {
//  double real_start_time = GetRealTime();
//  int n = 0, n1 = 0;
//  uint32_t i;
//  uint64_t sum = 0, len = 0;
//  fprintf(stderr, "[M::%s] kmer size: %d; skip: %d; #seq: %d\n", __func__, kmer_size_, window_size_, num_sequences);
//  for (i = 0; i < num_sequences; ++i) {
//    len += reference.GetSequenceLengthAt(i);
//  }
//  assert(len == reference.GetNumBases());
//  if (lookup_table_) {
//    n += kh_size(lookup_table_);
//  }
//  for (khint_t k = 0; k < kh_end(lookup_table_); ++k) {
//    if (kh_exist(lookup_table_, k)) {
//      sum += kh_key(lookup_table_, k) & 1 ? 1 : (uint32_t)kh_val(lookup_table_, k);
//      if (kh_key(lookup_table_, k) & 1) 
//        ++n1;
//    }
//  }
//  fprintf(stderr, "[M::%s::%.3f] distinct minimizers: %d (%.2f%% are singletons); average occurrences: %.3lf; average spacing: %.3lf\n",
//      __func__, GetRealTime() - real_start_time, n, 100.0*n1/n, (double)sum / n, (double)len / sum);
//}
//
//// always reserve space for minimizers in other functions
//void Index::GenerateMinimizerSketch(const SequenceBatch &sequence_batch, uint32_t sequence_index, std::vector<std::pair<uint64_t, uint64_t> > *minimizers) {
//  uint64_t num_shifted_bits = 2 * (kmer_size_ - 1); 
//  uint64_t mask = (((uint64_t)1) << (2 * kmer_size_)) - 1;
//  uint64_t seeds_in_two_strands[2] = {0, 0};
//  std::pair<uint64_t, uint64_t> buffer[256];
//  std::pair<uint64_t, uint64_t> min_seed = {UINT64_MAX, UINT64_MAX};
//  uint32_t sequence_length = sequence_batch.GetSequenceLengthAt(sequence_index);
//  const char *sequence = sequence_batch.GetSequenceAt(sequence_index);
//  assert(sequence_length > 0 && (window_size_ > 0 && window_size_ < 256) && (kmer_size_ > 0 && kmer_size_ <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
//  memset(buffer, 0xff, window_size_ * 16); // 2 uint64_t cost 16 bytes
//  int unambiguous_length = 0;
//  int position_in_buffer = 0;
//  int min_position = 0;
//  for (uint32_t position = 0; position < sequence_length; ++position) {
//    uint8_t current_base = SequenceBatch::CharToUint8(sequence[position]);
//    std::pair<uint64_t, uint64_t> current_seed = {UINT64_MAX, UINT64_MAX};
//    if (current_base < 4) { // not an ambiguous base
//      seeds_in_two_strands[0] = ((seeds_in_two_strands[0] << 2) | current_base) & mask; // forward k-mer
//      seeds_in_two_strands[1] = (seeds_in_two_strands[1] >> 2) | (((uint64_t)(3 ^ current_base)) << num_shifted_bits); // reverse k-mer
//      if (seeds_in_two_strands[0] == seeds_in_two_strands[1]) {
//        continue; // skip "symmetric k-mers" as we don't know it strand
//      }
//      uint64_t hash_keys_for_two_seeds[2] = {Hash64(seeds_in_two_strands[0], mask), Hash64(seeds_in_two_strands[1], mask)};
//      uint64_t strand = hash_keys_for_two_seeds[0] < hash_keys_for_two_seeds[1] ? 0 : 1; // strand
//      //uint64_t strand = seeds_in_two_strands[0] < seeds_in_two_strands[1] ? 0 : 1; // strand
//      ++unambiguous_length;
//      if (unambiguous_length >= kmer_size_) {
//        //current_seed.first = Hash64(seeds_in_two_strands[strand], mask);
//        current_seed.first = Hash64(hash_keys_for_two_seeds[strand], mask);
//        current_seed.second = ((((uint64_t)sequence_index) << 32 | (uint32_t)position) << 1) | strand;
//      }
//    } else {
//      unambiguous_length = 0;
//    }
//    buffer[position_in_buffer] = current_seed; // need to do this here as appropriate position_in_buffer and buf[position_in_buffer] are needed below
//    if (unambiguous_length == window_size_ + kmer_size_ - 1 && min_seed.first != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
//      for (int j = position_in_buffer + 1; j < window_size_; ++j)
//        if (min_seed.first == buffer[j].first && buffer[j].second != min_seed.second) 
//          minimizers->push_back(buffer[j]);
//      for (int j = 0; j < position_in_buffer; ++j)
//        if (min_seed.first == buffer[j].first && buffer[j].second != min_seed.second) 
//          minimizers->push_back(buffer[j]);
//    }
//    if (current_seed.first <= min_seed.first) { // a new minimum; then write the old min
//      if (unambiguous_length >= window_size_ + kmer_size_ && min_seed.first != UINT64_MAX) {
//        minimizers->push_back(min_seed);
//      }
//      min_seed = current_seed;
//      min_position = position_in_buffer;
//    } else if (position_in_buffer == min_position) { // old min has moved outside the window
//      if (unambiguous_length >= window_size_ + kmer_size_ - 1 && min_seed.first != UINT64_MAX) {
//        minimizers->push_back(min_seed);
//      }
//      min_seed.first = UINT64_MAX;
//      for (int j = position_in_buffer + 1; j < window_size_; ++j) { // the two loops are necessary when there are identical k-mers
//        if (min_seed.first >= buffer[j].first) {// >= is important s.t. min is always the closest k-mer
//          min_seed = buffer[j];
//          min_position = j; 
//        }
//      }
//      for (int j = 0; j <= position_in_buffer; ++j) {
//        if (min_seed.first >= buffer[j].first) {
//          min_seed = buffer[j];
//          min_position = j;
//        }
//      }
//      if (unambiguous_length >= window_size_ + kmer_size_ - 1 && min_seed.first != UINT64_MAX) { // write identical k-mers
//        for (int j = position_in_buffer + 1; j < window_size_; ++j) // these two loops make sure the output is sorted
//          if (min_seed.first == buffer[j].first && min_seed.second != buffer[j].second) 
//            minimizers->push_back(buffer[j]);
//        for (int j = 0; j <= position_in_buffer; ++j)
//          if (min_seed.first == buffer[j].first && min_seed.second != buffer[j].second) 
//            minimizers->push_back(buffer[j]);
//      }
//    }
//    ++position_in_buffer;
//    if (position_in_buffer == window_size_) {
//      position_in_buffer = 0;
//    }
//  }
//  if (min_seed.first != UINT64_MAX) {
//    minimizers->push_back(min_seed);
//  }
//}
//
//void Index::CheckIndex(uint32_t num_sequences, const SequenceBatch &reference) {
//  std::vector< std::pair<uint64_t, uint64_t> > tmp_table;
//  tmp_table.reserve(reference.GetNumBases() / window_size_ * 2);
//  for (uint32_t sequence_index = 0; sequence_index < num_sequences; ++sequence_index) {
//    GenerateMinimizerSketch(reference, sequence_index, &tmp_table);
//  }
//  std::cerr << "Collected " << tmp_table.size() << " minimizers.\n";
//  std::stable_sort(tmp_table.begin(), tmp_table.end());
//  std::cerr << "Sorted minimizers.\n";
//  uint32_t count = 0;
//  for (uint32_t i = 0; i < tmp_table.size(); ++i) {
//    khiter_t khash_iterator = kh_get(k64, lookup_table_, tmp_table[i].first << 1);
//    assert(khash_iterator != kh_end(lookup_table_));
//    uint64_t key = kh_key(lookup_table_, khash_iterator);
//    uint64_t value = kh_value(lookup_table_, khash_iterator);
//    if (key & 1) { //singleton
//      assert(tmp_table[i].second == value);
//      count = 0;
//    } else {
//      uint32_t offset = value >> 32;
//      uint32_t num_occ = value;
//      uint64_t value_in_index = occurrence_table_[offset + count];
//      assert(value_in_index == tmp_table[i].second);
//      ++count;
//      if (count == num_occ) {
//        count = 0;
//      }
//    }
//  }
//}

void SpatialIndex::GeneratePointCloudOnOneDirection(Direction direction, uint32_t signal_index, const std::vector<std::vector<bool> > &is_masked, const float *signal_values, size_t signal_length, int step_size, std::vector<Point> &point_cloud) {
  for (size_t signal_position = 0; signal_position < signal_length - dimension_ + 1; signal_position += step_size) {
    if (!is_masked[signal_index][signal_position]) {
      if (signal_position == 0 || (signal_position > 0 && abs(signal_values[signal_position] - point_cloud.back().value) > 0.01)) {
      uint64_t strand = direction == Positive ? 0 : 1;
      //uint64_t position = ((((uint64_t)signal_index) << 32 | (uint32_t)signal_position) << 1) | strand;
      uint64_t position = ((((uint64_t)signal_index) << 32 | (uint32_t)signal_position) << 1) | strand;
      point_cloud.emplace_back(position, signal_values[signal_position]);
      }
    }
  }
}

//void SpatialIndex::GetSignalIndexAndPosition(size_t point_index, size_t num_signals, const std::vector<std::vector<float> > &signals, size_t &signal_index, size_t &signal_position) {
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

void SpatialIndex::Construct(size_t num_signals, const std::vector<std::vector<bool> > &positive_is_masked, const std::vector<std::vector<bool> > &negative_is_masked, const std::vector<std::vector<float> > &positive_signals, const std::vector<std::vector<float> > &negative_signals) {
  double real_start_time = GetRealTime();
  int signal_point_step_size = 1;
  point_cloud_.reserve(2 * GetSignalsTotalLength(positive_signals));
  for (size_t signal_index = 0; signal_index < num_signals; ++signal_index) {
    GeneratePointCloudOnOneDirection(Positive, signal_index, positive_is_masked, positive_signals[signal_index].data(), positive_signals[signal_index].size(), signal_point_step_size, point_cloud_);
  }
  for (size_t signal_index = 0; signal_index < num_signals; ++signal_index) {
    GeneratePointCloudOnOneDirection(Negative, signal_index, negative_is_masked, negative_signals[signal_index].data(), negative_signals[signal_index].size(), signal_point_step_size, point_cloud_);
  }
  std::cerr << "Collected " << point_cloud_.size() << " points.\n";
  //sort(point_cloud_.begin(), point_cloud_.end());
  //std::cerr << "Sorted " << point_cloud_.size() << " points.\n";
  // TODO: mask repetitive kmers
  //std::cerr << "Masked repetitive kmers " << point_cloud_.size() << " points.\n";
  spatial_index_ = new SigmapAdaptor<float>(dimension_ /*dim*/, point_cloud_, max_leaf_ /* max leaf */);
  spatial_index_->index->buildIndex();
  std::cerr << "Built spatial index.\n";
  std::cerr << "Built index successfully in " << GetRealTime() - real_start_time << "s.\n";
}

void SpatialIndex::Save() {
  double real_start_time = GetRealTime();
  FILE *point_cloud_file = fopen((index_file_path_prefix_ + ".pt").c_str(), "wb");
  assert(point_cloud_file != NULL);
  fwrite(&dimension_, sizeof(int), 1, point_cloud_file);
  fwrite(&max_leaf_, sizeof(int), 1, point_cloud_file);
  size_t point_cloud_size = point_cloud_.size();
  fwrite(&point_cloud_size, sizeof(size_t), 1, point_cloud_file);
  // TODO(Haowen): we don't have to save the whole point cloud. Instead, we can just save the feature signal and build the point cloud on the fly. I will leave this as an easy later code refactor work.
  //for (size_t pi = 0; pi < point_cloud_size; ++pi) {
  //  fwrite(point_cloud_[pi].data(), sizeof(float), dimension_, point_cloud_file);
  //}
  fwrite(point_cloud_.data(), sizeof(Point), point_cloud_size, point_cloud_file);
  fclose(point_cloud_file);
  FILE *spatial_index_file = fopen((index_file_path_prefix_ + ".si").c_str(), "wb");
  assert(spatial_index_file != NULL);
  spatial_index_->index->saveIndex(spatial_index_file);
  fclose(spatial_index_file);
  std::cerr << "Saved in " << GetRealTime() - real_start_time << "s.\n";
}

void SpatialIndex::Load() {
  double real_start_time = GetRealTime();
  FILE *point_cloud_file = fopen((index_file_path_prefix_ + ".pt").c_str(), "rb");
  fread(&dimension_, sizeof(int), 1, point_cloud_file);
  fread(&max_leaf_, sizeof(int), 1, point_cloud_file);
  size_t point_cloud_size = 0;
  fread(&point_cloud_size, sizeof(size_t), 1, point_cloud_file);
  point_cloud_.resize(point_cloud_size);
  //for (size_t pi = 0; pi < point_cloud_size; ++pi) {
  //  point_cloud_[pi].resize(dimension_);
  //  fread(point_cloud_[pi].data(), sizeof(float), dimension_, point_cloud_file);
  //}
  fread(point_cloud_.data(), sizeof(Point), point_cloud_size, point_cloud_file);
  fclose(point_cloud_file);
  std::cerr << "Load point cloud successfully! dim: " << dimension_ << ", max leaf: " << max_leaf_ << ", point cloud size: " << point_cloud_size << "\n";
  FILE *spatial_index_file = fopen((index_file_path_prefix_ + ".si").c_str(), "rb");
  assert(spatial_index_file != NULL);
  spatial_index_ = new SigmapAdaptor<float>(dimension_ /*dim*/, point_cloud_, max_leaf_ /* max leaf */);
  spatial_index_->index->loadIndex(spatial_index_file);
  fclose(spatial_index_file);
  std::cerr << "Loaded index successfully in "<< GetRealTime() - real_start_time << "s.\n";
}

// TODO(Haowen): I need to modify the chaining algorithm of minimap2 to also consider the distance of points
//void SpatialIndex::GenerateCandidatesOnOneDirection(std::vector<uint64_t> *hits, std::vector<uint64_t> *candidates) {
//  hits->emplace_back(UINT64_MAX);
//  if (hits->size() > 0) {
//    std::sort(hits->begin(), hits->end());
//    int count = 1;
//    uint64_t previous_hit = (*hits)[0];
//    uint32_t previous_reference_id = previous_hit >> 32;
//    uint32_t previous_reference_position = previous_hit;
//    for (uint32_t pi = 1; pi < hits->size(); ++pi) {
//      uint32_t current_reference_id = (*hits)[pi] >> 32;
//      uint32_t current_reference_position = (*hits)[pi];
//      if (current_reference_id != previous_reference_id || current_reference_position - previous_reference_position > 1000) { // TODO(Haowen): find a proper parameter
//        if (count >= min_num_seeds_required_for_mapping_) {
//          candidates->push_back(previous_hit);
//        }
//        count = 1;
//      } else {
//        ++count;
//      }
//      previous_hit = (*hits)[pi];
//      previous_reference_id = current_reference_id;
//      previous_reference_position = current_reference_position;
//    }
//  }
//}
// //void SpatialIndex::CollectCandiates(int max_seed_frequency, const std::vector<std::vector<float> > &point_cloud, std::vector<uint64_t> *positive_hits, std::vector<uint64_t> *negative_hits) {
//  size_t num_points = point_cloud.size();
//  positive_hits->reserve(max_seed_frequencies_[0]);
//  negative_hits->reserve(max_seed_frequencies_[0]);
//  for (size_t pi = 0; pi < num_points; ++pi) {
//    size_t num_results = 5;
//    std::vector<size_t> ret_indexes(num_results);
//    std::vector<float> out_dists_sqr(num_results);
//    nanoflann::KNNResultSet<float> result_set(num_results);
//    result_set.init(&ret_indexes[0], &out_dists_sqr[0]);
//    spatial_index_->index->findNeighbors(result_set, point_cloud[pi].data(), nanoflann::SearchParams(10));
//    //std::cout << "knnSearch(nn=" <<num_results<< "): \n";
//    //for (size_t i = 0; i < num_results; i++) {
//    //  //size_t signal_index = 0;
//    //  //size_t signal_position = 0;
//    //  //GetSignalIndexAndPosition(ret_indexes[i], );
//    //  //std::cout << "ret_index["<<i<<"]=" << ret_indexes[i] << ", signal_index: " << ", signal_position: " << << " out_dist_sqr=" << out_dists_sqr[i] << std::endl;
//    //  std::cout << "ret_index["<<i<<"]=" << ret_indexes[i] << " out_dist_sqr=" << out_dists_sqr[i] << std::endl;
//    //}
//  } 
//  //for (uint32_t mi = 0; mi < num_minimizers; ++mi) {
//  //  khiter_t khash_iterator = kh_get(k64, lookup_table_, minimizers[mi].first << 1);
//  //  if (khash_iterator == kh_end(lookup_table_)) {
//  //    //std::cerr << "The minimizer is not in reference!\n";
//  //    continue;
//  //  }
//  //  uint64_t value = kh_value(lookup_table_, khash_iterator);
//  //  uint32_t read_position = minimizers[mi].second >> 1;
//  //  if (kh_key(lookup_table_, khash_iterator) & 1) { // singleton
//  //    uint64_t reference_id = value >> 33;
//  //    uint32_t reference_position = value >> 1;
//  //    // Check whether the strands of reference minimizer and read minimizer are the same
//  //    // Later, we can play some tricks with 0,1 here to make it faster.
//  //    if (((minimizers[mi].second & 1) ^ (value & 1)) == 0) { // same
//  //      uint32_t candidate_position = reference_position - read_position;// > 0 ? reference_position - read_position : 0;
//  //      // ok, for now we can't see the reference here. So let us don't do the check.
//  //      // Instead, we do it later some time when we check the candidates.
//  //      uint64_t candidate = (reference_id << 32) | candidate_position;
//  //      positive_hits->push_back(candidate);
//  //    } else {
//  //      uint32_t candidate_position = reference_position + read_position - kmer_size_ + 1;// < reference_length ? reference_position - read_position : 0;
//  //      uint64_t candidate = (reference_id << 32) | candidate_position;
//  //      negative_hits->push_back(candidate);
//  //    }
//  //  } else {
//  //    uint32_t offset = value >> 32;
//  //    uint32_t num_occurrences = value;
//  //    if (num_occurrences < (uint32_t)max_seed_frequency) {
//  //      for (uint32_t oi = 0; oi < num_occurrences; ++oi) {
//  //        uint64_t value = occurrence_table_[offset + oi];
//  //        uint64_t reference_id = value >> 33;
//  //        uint32_t reference_position = value >> 1;
//  //        if (((minimizers[mi].second & 1) ^ (value & 1)) == 0) { // same
//  //          uint32_t candidate_position = reference_position - read_position;
//  //          uint64_t candidate = (reference_id << 32) | candidate_position;
//  //          positive_hits->push_back(candidate);
//  //        } else {
//  //          uint32_t candidate_position = reference_position + read_position - kmer_size_ + 1;
//  //          uint64_t candidate = (reference_id << 32) | candidate_position;
//  //          negative_hits->push_back(candidate);
//  //        }
//  //      } 
//  //    }
//  //  }
//  //}
//}
//
//void SpatialIndex::GenerateCandidates(const std::vector<std::vector<float> > &point_cloud, std::vector<uint64_t> *positive_hits, std::vector<uint64_t> *negative_hits, std::vector<uint64_t> *positive_candidates, std::vector<uint64_t> *negative_candidates) {
//  CollectCandiates(max_seed_frequencies_[0], point_cloud, positive_hits, negative_hits);
//  // Now I can generate primer chain in candidates
//  // Let me use sort for now, but I can use merge later.
//  GenerateCandidatesOnOneDirection(positive_hits, positive_candidates);
//  GenerateCandidatesOnOneDirection(negative_hits, negative_candidates);
//  if (positive_candidates->size() + negative_candidates->size() == 0) {
//    positive_hits->clear();
//    negative_hits->clear();
//    CollectCandiates(max_seed_frequencies_[1], point_cloud, positive_hits, negative_hits);
//    GenerateCandidatesOnOneDirection(positive_hits, positive_candidates);
//    GenerateCandidatesOnOneDirection(negative_hits, negative_candidates);
//    // TODO: if necessary, we can further improve the rescue. But the code below is not thread safe. We can think about this later
////    if (positive_candidates->size() + negative_candidates->size() == 0) {
////      --min_num_seeds_required_for_mapping_;
////      min_num_seeds_required_for_mapping_ = std::max(min_num_seeds_required_for_mapping_, 1);
////      GenerateCandidatesOnOneDirection(positive_hits, positive_candidates);
////      GenerateCandidatesOnOneDirection(negative_hits, negative_candidates);
////    }
//  }
//}

void SpatialIndex::TacebackChains(int min_num_anchors, Direction direction, size_t chain_end_anchor_index, uint32_t chain_target_signal_index, const std::vector<float> &chaining_scores, const std::vector<size_t> &chaining_predecessors, const std::vector<std::vector<SignalAnchor> > &anchors_on_diff_signals, std::vector<bool> &anchor_is_used, std::vector<SignalAnchorChain> &chains) {
  if (!anchor_is_used[chain_end_anchor_index]) {
    //anchor_is_used[chain_end_anchor_index] = true;
    std::vector<SignalAnchor> anchors;
    anchors.reserve(100);
    bool stop_at_an_used_anchor = false;
    size_t chain_start_anchor_index = chain_end_anchor_index;//chaining_predecessors[chain_end_anchor_index];
    size_t chain_num_anchors = 1;
    //if (!anchor_is_used[chain_start_anchor_index] && chain_start_anchor_index != chain_end_anchor_index) {
    if (anchor_is_used[chaining_predecessors[chain_start_anchor_index]]) {
      stop_at_an_used_anchor = true;
    }
    while (!anchor_is_used[chaining_predecessors[chain_start_anchor_index]] && chaining_predecessors[chain_start_anchor_index] != chain_start_anchor_index) {
      anchors.push_back(anchors_on_diff_signals[chain_target_signal_index][chain_start_anchor_index]);
      anchor_is_used[chain_start_anchor_index] = true;
      chain_start_anchor_index = chaining_predecessors[chain_start_anchor_index];
      if (anchor_is_used[chaining_predecessors[chain_start_anchor_index]]) {
        stop_at_an_used_anchor = true;
      }
      ++chain_num_anchors;
    }
    //}
    if (chain_num_anchors >= (uint32_t)min_num_anchors) {
      float adjusted_chaining_score = chaining_scores[chain_end_anchor_index];
      if (stop_at_an_used_anchor) {
        adjusted_chaining_score -= chaining_scores[chain_start_anchor_index];
      }
      chains.emplace_back(SignalAnchorChain{adjusted_chaining_score, chain_target_signal_index, anchors_on_diff_signals[chain_target_signal_index][chain_start_anchor_index].target_position, anchors_on_diff_signals[chain_target_signal_index][chain_end_anchor_index].target_position, chain_num_anchors, 0, direction, anchors});
    }
  }
}

void SpatialIndex::GeneratePrimaryChains(std::vector<SignalAnchorChain> &chains) {
  std::sort(chains.begin(), chains.end(), std::greater<SignalAnchorChain>());
  //for (size_t i = 0; i < chains.size(); ++i) {
  //}
}

void SpatialIndex::ComputeMAPQ(std::vector<SignalAnchorChain> &chains) {
  if (chains.size() == 1) {
    chains[0].mapq = 60;
    return;
  } else {
    int mapq = 40 * (1 - chains[1].score / chains[0].score);// * std::min((size_t)1, chains[0].num_anchors / 20) * log(chains[0].score);
    if (mapq > 60) {
      mapq = 60;
    }
    if (mapq < 0) {
      mapq = 0;
    }
    chains[0].mapq = (uint8_t) mapq;
  }
}

void SpatialIndex::GenerateChains(const std::vector<float> &query_signal, int query_point_cloud_step_size, float search_radius, size_t num_target_signals, std::vector<SignalAnchorChain> &chains) {
  int max_gap_length = 5000; // TODO(Haowen): make it a parameter
  int chaining_band_length = 100; // TODO(Haowen): make it a parameter
  int min_num_anchors = 3; // TODO(Haowen): make it a parameter
  int num_best_chains = 2;
  int num_nearest_points = 500;
  std::vector<std::vector<std::vector<SignalAnchor> > > anchors_on_diff_signals(2);
  anchors_on_diff_signals[0] = std::vector<std::vector<SignalAnchor> >(num_target_signals);
  anchors_on_diff_signals[1] = std::vector<std::vector<SignalAnchor> >(num_target_signals);
  std::vector<std::vector<SignalAnchor> > &positive_anchors_on_diff_signals = anchors_on_diff_signals[0];
  std::vector<std::vector<SignalAnchor> > &negative_anchors_on_diff_signals = anchors_on_diff_signals[1];
  nanoflann::SearchParams params;
  params.sorted = true;
  std::vector<std::pair<size_t, float> > point_anchors;
  // Find reliable seeds
  std::vector<std::pair<float, size_t> > mean_diff_position;
  for (uint32_t pi = 0; pi < query_signal.size() - dimension_ + 1; ++pi) {
    float min_diff = std::numeric_limits<float>::max();
    for (int di = 1; di < dimension_; ++di) {
      float diff = abs(query_signal[pi + di] - query_signal[pi + di - 1]);
      if (diff < min_diff) {
        min_diff = diff;
      }
    }
    mean_diff_position.emplace_back(min_diff, pi);
  }
  std::sort(mean_diff_position.begin(), mean_diff_position.end(), compare1);
  // Collect anchors
  //for (uint32_t pi = 0; pi < query_signal.size() - dimension_ + 1; pi += query_point_cloud_step_size) {
  uint32_t previous_position = 0;
  uint32_t num_positions = 0;
  for (uint32_t pi = 0; pi < query_signal.size() - dimension_ + 1; ++pi) {
    uint32_t position = mean_diff_position[pi].second;
    //if (position <= previous_position + dimension_ / 2 && position + dimension_ / 2 >= previous_position) {
    //  continue;
    //}
    //std::cout << position << "\n";
    //size_t num_results = 100;
    //std::vector<size_t> ret_indexes(num_results);
    //std::vector<float> out_dists_sqr(num_results);
    //nanoflann::KNNResultSet<float> resultSet(num_results);
    //resultSet.init(&ret_indexes[0], &out_dists_sqr[0] );
    //size_t num_point_anchors = spatial_index_->index->findNeighbors(resultSet, query_signal.data() + position, nanoflann::SearchParams(10));
    size_t num_point_anchors = spatial_index_->index->radiusSearch(query_signal.data() + position, search_radius, point_anchors, params);
    //std::cout << "radiusSearch(): radius=" << search_radius << " -> " << num_point_anchors << " matches\n";
    //for (size_t ai = 0; ai < num_results; ai++) {
    float previous_distance = dimension_;
    for (size_t ai = 0; ai < num_point_anchors && ai < (uint32_t)num_nearest_points; ai++) {
      //if (point_anchors[ai].second > previous_distance * 2) {
      //  break;
      //} 
      //std::cout << "cloud size =" << point_cloud_.size() << " -> " << point_anchors[ai].first << " matches\n";
      Point &point = point_cloud_[point_anchors[ai].first];
      //Point &point = point_cloud_[ret_indexes[ai]];
      uint32_t target_signal_index = point.position >> 33, target_signal_position = point.position >> 1;
      Direction target_signal_direction = (point.position & 1) == 0 ? Positive : Negative;
      //GetSignalIndexAndPosition(point_anchors[ai].first, num_target_signals, target_signals, target_signal_index, target_signal_position);
      if (target_signal_direction == Positive) {
        positive_anchors_on_diff_signals[target_signal_index].emplace_back(SignalAnchor{target_signal_position, position, point_anchors[ai].second});
        //positive_anchors_on_diff_signals[target_signal_index].emplace_back(SignalAnchor{target_signal_position, position, out_dists_sqr[ai]});
      } else {
        negative_anchors_on_diff_signals[target_signal_index].emplace_back(SignalAnchor{target_signal_position, position, point_anchors[ai].second});
        //negative_anchors_on_diff_signals[target_signal_index].emplace_back(SignalAnchor{target_signal_position, position, out_dists_sqr[ai]});
      }
      previous_distance = point_anchors[ai].second;
      //std::cout << "idx["<< ai << "]=" << point_anchors[ai].first << " dist["<< ai << "]=" << point_anchors[ai].second << std::endl;
    }
    ++num_positions;
    if (num_positions >= (query_signal.size() - dimension_ + 1) / query_point_cloud_step_size) {
      break;
    }
    previous_position = position;
  }
  // Sort the anchors based on their occurrence on target signal
  for (size_t target_signal_index = 0; target_signal_index < num_target_signals; ++target_signal_index) {
    std::sort(positive_anchors_on_diff_signals[target_signal_index].begin(), positive_anchors_on_diff_signals[target_signal_index].end());
    std::sort(negative_anchors_on_diff_signals[target_signal_index].begin(), negative_anchors_on_diff_signals[target_signal_index].end());
  }
  // Chaining DP done on each individual target signal
  float max_chaining_score = std::numeric_limits<float>::min();
  //uint32_t max_chain_end_anchor_index = 0;
  //uint32_t max_chain_target_signal_index = 0;
  //uint32_t max_chain_num_anchors = 0;
  //uint32_t max_chain_start_anchor_index = 0;
  //Direction max_chain_direction = Positive;
  //bool max_chain_is_updated = false;
  for (size_t target_signal_index = 0; target_signal_index < num_target_signals; ++target_signal_index) {
    for (int direction_i = 0; direction_i < 2; ++direction_i) {
      //max_chain_is_updated = false;
      std::vector<float> chaining_scores;
      std::vector<size_t> chaining_predecessors;
      std::vector<bool> anchor_is_used;
      std::vector<std::pair<float, size_t> > end_anchor_index_chaining_scores;
      for (size_t anchor_index = 0; anchor_index < anchors_on_diff_signals[direction_i][target_signal_index].size(); ++anchor_index) {
        float distance_coefficient = 1 - 0.2 * anchors_on_diff_signals[direction_i][target_signal_index][anchor_index].distance / search_radius;
        chaining_scores.emplace_back(distance_coefficient * dimension_);
        chaining_predecessors.emplace_back(anchor_index);
        anchor_is_used.emplace_back(false);
        int32_t current_anchor_target_position = anchors_on_diff_signals[direction_i][target_signal_index][anchor_index].target_position;
        int32_t current_anchor_query_position = anchors_on_diff_signals[direction_i][target_signal_index][anchor_index].query_position;
        size_t previous_anchor_index = 0;
        if (anchor_index > (size_t)chaining_band_length) {
          previous_anchor_index = anchor_index - chaining_band_length;
        }
        for (; previous_anchor_index < anchor_index; ++previous_anchor_index) {
          int32_t previous_anchor_target_position = anchors_on_diff_signals[direction_i][target_signal_index][previous_anchor_index].target_position;
          int32_t previous_anchor_query_position = anchors_on_diff_signals[direction_i][target_signal_index][previous_anchor_index].query_position;
          if (previous_anchor_query_position == current_anchor_query_position) {
            continue;
          }
          if (previous_anchor_target_position == current_anchor_target_position) {
            continue;
          }
          int32_t target_position_diff = current_anchor_target_position - previous_anchor_target_position;
          assert(target_position_diff > 0);
          int32_t query_position_diff = current_anchor_query_position - previous_anchor_query_position;
          float current_chaining_score = 0;
          //std::cerr << "curr_ai: " << anchor_index << ", adjusted_d: " << dimension_ * distance_coefficient << ", pre_ai: " << previous_anchor_index << ", ";
          //std::cerr << "target_diff: " << target_position_diff << ", query_diff: " << query_position_diff << "\n";
          if (query_position_diff < 0) {
            current_chaining_score = std::numeric_limits<float>::min();
          } else {
            float matching_dimensions = std::min(std::min(target_position_diff, query_position_diff), dimension_) * distance_coefficient;
            int gap_length = std::abs(target_position_diff - query_position_diff);
            float gap_scale = target_position_diff > 0 ?(float)query_position_diff / target_position_diff : 1; 
            float gap_cost = 0;
            if (gap_length != 0) {
              //if (gap_length < max_gap_length && target_position_diff <= 1000 && query_position_diff <= 1000 * gap_scale && gap_scale < 1.8 && gap_scale > 0.95) { 
              if (gap_length < max_gap_length && gap_scale < 1.8 && gap_scale > 0.90) { 
                //if (gap_length > target_position_diff) {
                  //gap_cost = 0.001 * dimension_ * distance_coefficient * gap_length + 0.05 * std::log2(gap_length);
                //if (gap_scale > 1) {
                //gap_cost = exp(-1.2 / gap_scale);
                //}
                //}
                current_chaining_score = chaining_scores[previous_anchor_index] + matching_dimensions - gap_cost;
              } else {
                gap_cost = std::numeric_limits<float>::max();
                current_chaining_score = std::numeric_limits<float>::min();//chaining_scores[previous_anchor_index] + matching_dimensions - gap_cost;
              }
            }
            //if (gap_length != 0) {
            //  if (gap_length < max_gap_length && 0.9 * target_position_diff <= query_position_diff && query_position_diff < 2 * target_position_diff) { 
            //    if (gap_length > target_position_diff) {
            //      //gap_cost = 0.001 * dimension_ * distance_coefficient * gap_length + 0.05 * std::log2(gap_length);
            //      gap_cost = 0.25 * std::log2(gap_length);
            //    }
            //    current_chaining_score = chaining_scores[previous_anchor_index] + matching_dimensions - gap_cost;
            //  } else {
            //    gap_cost = std::numeric_limits<float>::max();
            //    current_chaining_score = std::numeric_limits<float>::min();//chaining_scores[previous_anchor_index] + matching_dimensions - gap_cost;
            //  }
            //}
            //std::cerr << ", matching_d: " << matching_dimensions;
            //std::cerr << ", gap_len: " << gap_length;
            //std::cerr << ", gap_cost: " << gap_cost;
            //std::cerr << ", current chaining score: " << current_chaining_score << "\n";
          }
          if (current_chaining_score > chaining_scores[anchor_index]) {
            chaining_scores[anchor_index] = current_chaining_score;
            chaining_predecessors[anchor_index] = previous_anchor_index;
            //std::cerr << "update_curr_max: " << current_chaining_score << "\n";
          }
        }
        // Update chain with max score
        if (chaining_scores[anchor_index] > max_chaining_score) {
          max_chaining_score = chaining_scores[anchor_index];
        //  max_chain_end_anchor_index = anchor_index;
        //  max_chain_target_signal_index = target_signal_index;
        //  max_chain_direction = direction_i == 0 ? Positive : Negative;
        //  max_chain_is_updated = true;
        }
        if (chaining_scores.back() >= max_chaining_score / 2) {
          end_anchor_index_chaining_scores.emplace_back(chaining_scores.back(), anchor_index);
        }
      }
      // Sort a vector of <end anchor index, chaining score> 
      std::sort(end_anchor_index_chaining_scores.begin(), end_anchor_index_chaining_scores.end(), compare);
      // Traceback all chains from higest score to lowest
      for (size_t anchor_index = 0; anchor_index < end_anchor_index_chaining_scores.size() && anchor_index < (size_t)num_best_chains; ++anchor_index) {
        TacebackChains(min_num_anchors, direction_i == 0 ? Positive : Negative, end_anchor_index_chaining_scores[anchor_index].second, target_signal_index, chaining_scores, chaining_predecessors, anchors_on_diff_signals[direction_i], anchor_is_used, chains);
        if (chaining_scores[end_anchor_index_chaining_scores[anchor_index].second] < max_chaining_score / 2) {
          break;
        }
      }
      //if (max_chain_is_updated) {
      //  max_chain_start_anchor_index = chaining_predecessors[max_chain_end_anchor_index];
      //  if (max_chain_start_anchor_index != chaining_predecessors[max_chain_start_anchor_index]) {
      //    max_chain_num_anchors = 1;
      //    while (chaining_predecessors[max_chain_start_anchor_index] != max_chain_start_anchor_index) {
      //      max_chain_start_anchor_index = chaining_predecessors[max_chain_start_anchor_index];
      //      ++max_chain_num_anchors;
      //    }
      //  }
      //}
    }
  }
  if (chains.size() > 0) {
    // Generate primary chains
    GeneratePrimaryChains(chains);
    // Compute MAPQ
    ComputeMAPQ(chains);
  }
  //if (max_chain_direction == Positive) {
  //  positive_chains.emplace_back(SignalAnchorChain{max_chaining_score, max_chain_target_signal_index, positive_anchors_on_diff_signals[max_chain_target_signal_index][max_chain_start_anchor_index].target_position, positive_anchors_on_diff_signals[max_chain_target_signal_index][max_chain_end_anchor_index].target_position, max_chain_num_anchors});
  //} else {
  //  negative_chains.emplace_back(SignalAnchorChain{max_chaining_score, max_chain_target_signal_index, negative_anchors_on_diff_signals[max_chain_target_signal_index][max_chain_start_anchor_index].target_position, negative_anchors_on_diff_signals[max_chain_target_signal_index][max_chain_end_anchor_index].target_position, max_chain_num_anchors});
  //}
  //std::cerr << "Max chaining score: " << max_chaining_score << ", signal_index: " << max_chain_target_signal_index << ", anchor target end postion: " << anchors_on_diff_signals[max_chain_target_signal_index].target_position << ".\n"
}
} // namespace sigmap
