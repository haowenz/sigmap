#include "sigmap.h"

#include <cassert>
#include <iostream>
#include <omp.h>
#include <string>

#include "cxxopts.hpp"
#include "khash.h"
#include "nanoflann.hpp"
#include "spatial_index.h"
#include "sequence_batch.h"
#include "pore_model.h"

KHASH_MAP_INIT_INT(k64, uint64_t);

namespace sigmap {
void Sigmap::GenerateMaskedPositions(int kmer_size, float frequency, uint32_t num_sequences, const SequenceBatch &sequence_batch, std::vector<std::vector<bool> > &positive_is_masked, std::vector<std::vector<bool> > &negative_is_masked) {
  double real_start_time = GetRealTime();
  khash_t(k64)* kmer_hist = kh_init(k64);
  uint64_t num_shifted_bits = 2 * (kmer_size - 1); 
  uint64_t mask = (((uint64_t)1) << (2 * kmer_size)) - 1;
  uint64_t seeds_in_two_strands[2] = {0, 0};
  uint64_t num_kmers = 0;
  int unambiguous_length = 0;
  // Count kmers
  for (uint32_t sequence_index = 0; sequence_index < num_sequences; ++sequence_index) {
    uint32_t sequence_length = sequence_batch.GetSequenceLengthAt(sequence_index);
    const char *sequence = sequence_batch.GetSequenceAt(sequence_index);
    unambiguous_length = 0;
    seeds_in_two_strands[0] = 0;
    seeds_in_two_strands[1] = 0;
    for (uint32_t position = 0; position < sequence_length; ++position) {
      uint8_t current_base = SequenceBatch::CharToUint8(sequence[position]);
      if (current_base < 4) { // not an ambiguous base
        seeds_in_two_strands[0] = ((seeds_in_two_strands[0] << 2) | current_base) & mask; // forward k-mer
        seeds_in_two_strands[1] = (seeds_in_two_strands[1] >> 2) | (((uint64_t)(3 ^ current_base)) << num_shifted_bits); // reverse k-mer
        //if (seeds_in_two_strands[0] == seeds_in_two_strands[1]) {
        //  continue; // skip "symmetric k-mers" as we don't know it strand
        //}
        ++unambiguous_length;
        if (unambiguous_length >= kmer_size) {
          uint64_t strand = seeds_in_two_strands[0] < seeds_in_two_strands[1] ? 0 : 1; // strand
          khiter_t kmer_hist_iterator = kh_get(k64, kmer_hist, seeds_in_two_strands[strand]);
          if (kmer_hist_iterator != kh_end(kmer_hist)) {
            kh_value(kmer_hist, kmer_hist_iterator) += 1;
          } else {
            int khash_return_code;
            khiter_t kmer_hist_insert_iterator = kh_put(k64, kmer_hist, seeds_in_two_strands[strand], &khash_return_code);
            assert(khash_return_code != -1 && khash_return_code != 0);
            kh_value(kmer_hist, kmer_hist_insert_iterator) = 1;
          }
          ++num_kmers;
        }
      } else {
        unambiguous_length = 0;
        seeds_in_two_strands[0] = 0;
        seeds_in_two_strands[1] = 0;
      }
    }
  }
  std::cout << "# " << kmer_size << "mers: " << num_kmers << "\n";
  // Mask kmers
  uint64_t num_masked_kmers = 0;
  for (uint32_t sequence_index = 0; sequence_index < num_sequences; ++sequence_index) {
    uint32_t sequence_length = sequence_batch.GetSequenceLengthAt(sequence_index);
    const char *sequence = sequence_batch.GetSequenceAt(sequence_index);
    positive_is_masked.emplace_back(std::vector<bool>(sequence_length - kmer_size + 1, false));
    unambiguous_length = 0;
    seeds_in_two_strands[0] = 0;
    seeds_in_two_strands[1] = 0;
    for (uint32_t position = 0; position < sequence_length; ++position) {
      uint8_t current_base = SequenceBatch::CharToUint8(sequence[position]);
      if (current_base < 4) { // not an ambiguous base
        seeds_in_two_strands[0] = ((seeds_in_two_strands[0] << 2) | current_base) & mask; // forward k-mer
        seeds_in_two_strands[1] = (seeds_in_two_strands[1] >> 2) | (((uint64_t)(3 ^ current_base)) << num_shifted_bits); // reverse k-mer
        ++unambiguous_length;
        if (unambiguous_length >= kmer_size) {
          uint64_t strand = seeds_in_two_strands[0] < seeds_in_two_strands[1] ? 0 : 1; // strand
          khiter_t kmer_hist_iterator = kh_get(k64, kmer_hist, seeds_in_two_strands[strand]);
          assert(kmer_hist_iterator != kh_end(kmer_hist));
          uint64_t kmer_freq = kh_value(kmer_hist, kmer_hist_iterator);
          //std::cerr << "count: " << kh_value(kmer_hist, kmer_hist_iterator) << " position: " << position << " sequence_i: " << sequence_index << "\n";
          //std::cerr << "size: " << is_masked[sequence_index].size() << "\n";
          if ((kmer_freq / (float)num_kmers) > frequency) {
            positive_is_masked[sequence_index][position + 1 - kmer_size] = true;
            ++num_masked_kmers;
          } else {
            positive_is_masked[sequence_index][position + 1 - kmer_size] = false;
          }
        }
      } else {
        unambiguous_length = 0;
        seeds_in_two_strands[0] = 0;
        seeds_in_two_strands[1] = 0;
        if (position >= (uint32_t)kmer_size - 1) {
          positive_is_masked[sequence_index][position + 1 - kmer_size] = true;
          ++num_masked_kmers;
        }
      }
    }
    const char *negative_sequence = sequence_batch.GetNegativeSequenceAt(sequence_index).data();
    negative_is_masked.emplace_back(std::vector<bool>(sequence_length - kmer_size + 1, false));
    unambiguous_length = 0;
    seeds_in_two_strands[0] = 0;
    seeds_in_two_strands[1] = 0;
    for (uint32_t position = 0; position < sequence_length; ++position) {
      uint8_t current_base = SequenceBatch::CharToUint8(negative_sequence[position]);
      if (current_base < 4) { // not an ambiguous base
        seeds_in_two_strands[0] = ((seeds_in_two_strands[0] << 2) | current_base) & mask; // forward k-mer
        seeds_in_two_strands[1] = (seeds_in_two_strands[1] >> 2) | (((uint64_t)(3 ^ current_base)) << num_shifted_bits); // reverse k-mer
        ++unambiguous_length;
        if (unambiguous_length >= kmer_size) {
          uint64_t strand = seeds_in_two_strands[0] < seeds_in_two_strands[1] ? 0 : 1; // strand
          khiter_t kmer_hist_iterator = kh_get(k64, kmer_hist, seeds_in_two_strands[strand]);
          assert(kmer_hist_iterator != kh_end(kmer_hist));
          uint64_t kmer_freq = kh_value(kmer_hist, kmer_hist_iterator);
          //std::cerr << "count: " << kh_value(kmer_hist, kmer_hist_iterator) << " position: " << position << " sequence_i: " << sequence_index << "\n";
          //std::cerr << "size: " << is_masked[sequence_index].size() << "\n";
          if ((kmer_freq / (float)num_kmers) > frequency) {
            negative_is_masked[sequence_index][position + 1 - kmer_size] = true;
            ++num_masked_kmers;
          } else {
            negative_is_masked[sequence_index][position + 1 - kmer_size] = false;
          }
        }
      } else {
        unambiguous_length = 0;
        seeds_in_two_strands[0] = 0;
        seeds_in_two_strands[1] = 0;
        if (position >= (uint32_t)kmer_size - 1) {
          negative_is_masked[sequence_index][position + 1 - kmer_size] = true;
          ++num_masked_kmers;
        }
      }
    }
  }
  kh_destroy(k64, kmer_hist);
  std::cerr << "Mask " << num_masked_kmers << " high frequency kmer in " << GetRealTime() - real_start_time << "s.\n";
}

//void Sigmap::EmplaceBackMappingRecord(uint32_t read_id, const char *read_name, uint32_t read_length, uint32_t read_start_position, uint32_t read_end_position, uint32_t barcode, uint32_t fragment_start_position, uint32_t fragment_length, uint8_t mapq, uint8_t direction, uint8_t is_unique, std::vector<PAFMapping> *mappings_on_diff_ref_seqs) {
//  mappings_on_diff_ref_seqs->emplace_back(PAFMapping{read_id, std::string(read_name), read_length, read_start_position, read_end_position, fragment_start_position, fragment_length, mapq, direction, is_unique});
//}

void Sigmap::OutputMappingsInVector(uint8_t mapq_threshold, uint32_t num_reference_sequences, const SequenceBatch &reference, const std::vector<std::vector<PAFMapping> > &mappings) {
  for (uint32_t ri = 0; ri < num_reference_sequences; ++ri) {
    for (auto it = mappings[ri].begin(); it != mappings[ri].end(); ++it) {
      uint8_t mapq = (it->mapq);
      //uint8_t is_unique = (it->is_unique);
      if (mapq >= mapq_threshold && mapq <= 60) {
        //if (allocate_multi_mappings_ || (only_output_unique_mappings_ && is_unique == 1)) {
          output_tools_->AppendMapping(ri, reference, *it);
        //}
      } else {
          output_tools_->AppendUnmappedRead(ri, reference, *it);
      }
    }
  }
}

uint32_t Sigmap::MoveMappingsInBuffersToMappingContainer(uint32_t num_reference_sequences, std::vector<std::vector<std::vector<PAFMapping> > > *mappings_on_diff_ref_seqs_for_diff_threads_for_saving) {
  double real_start_time = GetRealTime();
  uint32_t num_moved_mappings = 0;
  for (int ti = 0; ti < num_threads_; ++ti) {
    for (uint32_t i = 0; i < num_reference_sequences; ++i) {
      num_moved_mappings += (*mappings_on_diff_ref_seqs_for_diff_threads_for_saving)[ti][i].size();
      mappings_on_diff_ref_seqs_[i].insert(mappings_on_diff_ref_seqs_[i].end(), std::make_move_iterator((*mappings_on_diff_ref_seqs_for_diff_threads_for_saving)[ti][i].begin()), std::make_move_iterator((*mappings_on_diff_ref_seqs_for_diff_threads_for_saving)[ti][i].end()));
      (*mappings_on_diff_ref_seqs_for_diff_threads_for_saving)[ti][i].clear();
    }
  }
  std::cerr << "Move mappings in " << GetRealTime() - real_start_time << "s.\n";
  return num_moved_mappings;
}

void Sigmap::Map() {
  // Load read signals
  SignalBatch read_signal_batch;
  read_signal_batch.InitializeLoading(signal_directory_);
  size_t num_loaded_read_signals = read_signal_batch.LoadAllReadSignals();
  // Load pore model
  PoreModel pore_model;
  pore_model.Load(pore_model_file_path_);
  // Load reference genome
  SequenceBatch reference_sequence_batch;
  reference_sequence_batch.InitializeLoading(reference_file_path_);
  uint32_t num_reference_sequences = reference_sequence_batch.LoadAllSequences();
  // Get reverse complement of each ref seq
  for (size_t reference_sequence_index = 0; reference_sequence_index < num_reference_sequences; ++reference_sequence_index) {
    reference_sequence_batch.PrepareNegativeSequenceAt(reference_sequence_index);
  }
  // Use pore model to convert reference sequence to signal
  SignalBatch reference_signal_batch;
  reference_signal_batch.ConvertSequencesToSignals(reference_sequence_batch, pore_model, num_reference_sequences);
  // Normalize reference signals
  std::vector<std::vector<float> > positive_reference_feature_signals;
  std::vector<std::vector<float> > negative_reference_feature_signals;
  for (size_t reference_signal_index = 0; reference_signal_index < num_reference_sequences; ++reference_signal_index) {
    positive_reference_feature_signals.push_back(std::vector<float>());
    negative_reference_feature_signals.push_back(std::vector<float>());
    GenerateZscoreNormalizedSignal(reference_signal_batch.GetSignalAt(reference_signal_index).signal_values, reference_signal_batch.GetSignalLengthAt(reference_signal_index), positive_reference_feature_signals.back());
    GenerateZscoreNormalizedSignal(reference_signal_batch.GetSignalAt(reference_signal_index).negative_signal_values, reference_signal_batch.GetSignalLengthAt(reference_signal_index), negative_reference_feature_signals.back());
  }
  // Load spatial index for reference signals 
  SpatialIndex reference_spatial_index(1000, std::vector<int>(1000,5000), reference_index_file_path_);
  reference_spatial_index.Load();

  mappings_on_diff_ref_seqs_.reserve(num_reference_sequences);
  for (uint32_t i = 0; i < num_reference_sequences; ++i) {
    mappings_on_diff_ref_seqs_.emplace_back(std::vector<PAFMapping>());
  }
  output_tools_ = std::unique_ptr<PAFOutputTools<PAFMapping> >(new PAFOutputTools<PAFMapping>);
  //std::vector<std::vector<PAFMapping> > &mappings;
  std::vector<std::vector<std::vector<PAFMapping> > > mappings_on_diff_ref_seqs_for_diff_threads;
  mappings_on_diff_ref_seqs_for_diff_threads.reserve(num_threads_);
  //mappings_on_diff_ref_seqs_for_diff_threads_for_saving.reserve(num_threads_);
  for (int ti = 0; ti < num_threads_; ++ti) {
    mappings_on_diff_ref_seqs_for_diff_threads.emplace_back(std::vector<std::vector<PAFMapping> >(num_reference_sequences));
    //mappings_on_diff_ref_seqs_for_diff_threads_for_saving.emplace_back(std::vector<std::vector<MappingRecord> >(num_reference_sequences));
    for (uint32_t i = 0; i < num_reference_sequences; ++i) {
      mappings_on_diff_ref_seqs_for_diff_threads[ti][i].reserve(5000 / num_threads_ / num_reference_sequences);
      //mappings_on_diff_ref_seqs_for_diff_threads[ti][i].reserve((num_loaded_pairs + num_loaded_pairs / 1000 * max_num_best_mappings_) / num_threads_ / num_reference_sequences);
      //mappings_on_diff_ref_seqs_for_diff_threads_for_saving[ti][i].reserve((num_loaded_pairs + num_loaded_pairs / 1000 * max_num_best_mappings_) / num_threads_ / num_reference_sequences);
    }
  }
  output_tools_->InitializeMappingOutput(output_file_path_);
  //output_tools_->AppendMapping(last_rid, reference, last_mapping);

  // Map each reads
  double real_start_time = GetRealTime();
#pragma omp parallel default(none) shared(reference_sequence_batch, positive_reference_feature_signals, negative_reference_feature_signals, reference_spatial_index, read_signal_batch, std::cerr, num_loaded_read_signals, num_reference_sequences, mappings_on_diff_ref_seqs_for_diff_threads) num_threads(num_threads_)
  {
  std::vector<float> read_feature_signal;
  std::vector<Point> read_point_cloud;
  std::vector<SignalAnchorChain> chains;
#pragma omp single
  {
  //int grain_size = num_loaded_read_signals / num_threads_ / 100;
  //grain_size = grain_size > 0 ? grain_size : 1;
#pragma omp taskloop //grainsize(grain_size) //num_tasks(num_threads_* 50)
  for (size_t read_signal_index = 0; read_signal_index < num_loaded_read_signals; ++read_signal_index) {
    double real_mapping_start_time = GetRealTime();
    std::cerr << "mapping read " << read_signal_index << ".\n";
    read_feature_signal.clear();
    //read_signal_batch.MovingMedianSignalAt(read_signal_index, 8);
    //read_signal_batch.NormalizeSignalAt(read_signal_index);
    GenerateEvents(read_signal_batch.GetSignalAt(read_signal_index), read_feature_signal);
    //if (read_feature_signal.size() > 2000) {
    //  read_feature_signal.erase(read_feature_signal.begin() + 2000, read_feature_signal.end());
    //}
    if (read_feature_signal.size() > 50) {
      read_point_cloud.clear();
      chains.clear();
      float search_radius = 0.08;
      //float search_radius = 0.05;
      int read_signal_point_cloud_step_size = 8;
      //if (read_feature_signal.size() < 10000) {
      //  read_signal_point_cloud_step_size = 7;
      //} 
      //if (read_feature_signal.size() < 8000) {
      //  read_signal_point_cloud_step_size = 5;
      //}
      //if (read_feature_signal.size() < 5000) {
      //  read_signal_point_cloud_step_size = 3;
      //}
      //if (read_feature_signal.size() < 2000) {
      //  read_signal_point_cloud_step_size = 1;
      //}
      read_signal_point_cloud_step_size = 1;
      reference_spatial_index.GenerateChains(read_feature_signal, read_signal_point_cloud_step_size, search_radius, num_reference_sequences, chains);
      // Save results in vector and output PAF
      double mapping_time = GetRealTime() - real_mapping_start_time;
      std::vector<std::vector<PAFMapping> > &mappings_on_diff_ref_seqs = mappings_on_diff_ref_seqs_for_diff_threads[omp_get_thread_num()];
      if (chains.size() > 0) {
        float anchor_ref_gap_avg_length = 0;
        float anchor_read_gap_avg_length = 0;
        float average_anchor_distance = 0;
        for (size_t ai = 0; ai < chains[0].anchors.size(); ++ai) {
          average_anchor_distance += chains[0].anchors[ai].distance;
          if (ai < chains[0].anchors.size() - 1) {
            anchor_ref_gap_avg_length += chains[0].anchors[ai].target_position - chains[0].anchors[ai + 1].target_position;
            anchor_read_gap_avg_length += chains[0].anchors[ai].query_position - chains[0].anchors[ai + 1].query_position;
          }
          //std::cerr << "(" << chains[i].anchors[ai].target_position << "," << chains[i].anchors[ai].query_position << "," << chains[i].anchors[ai].distance << "), ";
        }
        average_anchor_distance /= chains[0].num_anchors;
        anchor_ref_gap_avg_length /= chains[0].num_anchors;
        anchor_read_gap_avg_length /= chains[0].num_anchors;
        std::string tags;
        tags.append("mt:f:" + std::to_string(mapping_time * 1000));
        tags.append("\tsl:i:" + std::to_string(read_signal_batch.GetSignalLengthAt(read_signal_index)));
        tags.append("\tcm:i:" + std::to_string(chains[0].num_anchors));
        tags.append("\ts1:f:" + std::to_string(chains[0].score));
        tags.append("\ts2:f:" + std::to_string(chains.size() > 1 ? chains[1].score : 0));
        tags.append("\tad:f:" + std::to_string(average_anchor_distance));
        tags.append("\tat:f:" + std::to_string(anchor_ref_gap_avg_length));
        tags.append("\taq:f:" + std::to_string(anchor_read_gap_avg_length));
        mappings_on_diff_ref_seqs[chains[0].reference_sequence_index].emplace_back(PAFMapping{(uint32_t)read_signal_index, std::string(read_signal_batch.GetSignalNameAt(read_signal_index)), (uint32_t)read_feature_signal.size(), chains[0].anchors.back().query_position, chains[0].anchors[0].query_position, chains[0].direction == Positive ? (uint32_t)chains[0].start_position : (uint32_t)(reference_sequence_batch.GetSequenceLengthAt(chains[0].reference_sequence_index) + 1 - chains[0].end_position), (uint32_t)(chains[0].end_position - chains[0].start_position + 1), chains[0].mapq, chains[0].direction == Positive ? (uint8_t)1 : (uint8_t)0, (uint8_t)1, tags});

        //EmplaceBackMappingRecord(read_signal_index, read_signal_batch.GetSignalNameAt(read_signal_index), read_signal_batch.GetSignalLengthAt(read_signal_index), chains[0].direction == Positive ? 0 : 1, chains[0].direction == Positive ? chains[0].start_position : reference_sequence_batch.GetSequenceLengthAt(chains[0].reference_sequence_index) + 1 - chains[0].end_position, chains[0].end_position - chains[0].start_position + 1, chains[0].mapq, chains[0].direction == Positive ? 1 : 0, 1, &(mappings_on_diff_ref_seqs[chains[0].reference_sequence_index]));
#ifdef DEBUG
        //if (chains[0].direction == Positive) {
        //  std::cerr << "Direction: positive.\n";
        //} else {
        //  std::cerr << "Direction: negative.\n";
        //}
        //std::cerr << "Best chaining score: " << chains[0].score << ", signal_index: " << chains[0].reference_sequence_index << ", anchor target start postion: " << chains[0].start_position << ", anchor target end postion: " << chains[0].end_position << ", # anchors: " << chains[0].num_anchors << ", mapq: " << (int)chains[0].mapq << ".\n";
        for (size_t i = 0; i < chains.size(); ++i) {
          if (chains[i].direction == Positive) {
            std::cerr << i << " best chaining score: " << chains[i].score << ", direction: +" << ", reference name: " << reference_sequence_batch.GetSequenceNameAt(chains[i].reference_sequence_index) << ", anchor target start postion: " << chains[i].start_position << ", anchor target end postion: " << chains[i].end_position << ", # anchors: " << chains[i].num_anchors << ", mapq: " << (int)chains[i].mapq << ".\n";
            for (size_t ai = 0; ai < chains[i].anchors.size(); ++ai) {
              std::cerr << "(" << chains[i].anchors[ai].target_position << "," << chains[i].anchors[ai].query_position << "," << chains[i].anchors[ai].distance << "), ";
            }
          } else {
            std::cerr << i << " best chaining score: " << chains[i].score << ", direction: -" << ", reference name: " << reference_sequence_batch.GetSequenceNameAt(chains[i].reference_sequence_index) << ", anchor target start postion: " << reference_sequence_batch.GetSequenceLengthAt(chains[i].reference_sequence_index) + 1 - chains[i].end_position << ", anchor target end postion: " << reference_sequence_batch.GetSequenceLengthAt(chains[i].reference_sequence_index) - chains[i].start_position << ", # anchors: " << chains[i].num_anchors << ", mapq: " << (int)chains[i].mapq << ".\n";
            for (size_t ai = 0; ai < chains[i].anchors.size(); ++ai) {
              std::cerr << "(" << chains[i].anchors[ai].target_position << "," << chains[i].anchors[ai].query_position << "," << chains[i].anchors[ai].distance << "), ";
            }
          }
          std::cerr << "\n";
          std::cerr << "\n";
        }
        std::cerr << "Read name: " << read_signal_batch.GetSignalNameAt(read_signal_index) << ", length: " << read_feature_signal.size() << ", reference name: " << reference_sequence_batch.GetSequenceNameAt(chains[0].reference_sequence_index) << ", length: " << positive_reference_feature_signals[chains[0].reference_sequence_index].size() << "\n";
        std::cerr << "\n";
#endif
      } else {
        std::string tags;
        tags.append("mt:f:" + std::to_string(mapping_time * 1000));
        tags.append("\tsl:i:" + std::to_string(read_signal_batch.GetSignalLengthAt(read_signal_index)));
        tags.append("\tcm:i:" + std::to_string(0));
        tags.append("\ts1:f:" + std::to_string(0));
        tags.append("\ts2:f:" + std::to_string(0));
        mappings_on_diff_ref_seqs[0].emplace_back(PAFMapping{(uint32_t)read_signal_index, std::string(read_signal_batch.GetSignalNameAt(read_signal_index)), (uint32_t)read_feature_signal.size(), 0, 0, 0, 0, 61, (uint8_t)0, (uint8_t)1, tags});
      }
    }
  }
  } // end of openmp single
  } // end of openmp parallel
  std::cerr << "Finished mapping in " << GetRealTime() - real_start_time << ", # reads: " << num_loaded_read_signals << "\n";

  MoveMappingsInBuffersToMappingContainer(num_reference_sequences, &mappings_on_diff_ref_seqs_for_diff_threads);
  OutputMappingsInVector(0, num_reference_sequences, reference_sequence_batch, mappings_on_diff_ref_seqs_);
  output_tools_->FinalizeMappingOutput();
  read_signal_batch.FinalizeLoading();
  reference_sequence_batch.FinalizeLoading();
  reference_signal_batch.FinalizeLoading();
}

void Sigmap::DTWAlign() {
  SignalBatch read_signal_batch;
  read_signal_batch.InitializeLoading(signal_directory_);
  size_t num_loaded_read_signals = read_signal_batch.LoadAllReadSignals();
  double real_normalization_start_time = GetRealTime();
  for (size_t read_index = 0; read_index < num_loaded_read_signals; ++read_index) {
  read_signal_batch.NormalizeSignalAt(read_index);
  }
  std::cerr << "Normalize " << num_loaded_read_signals << " read signals in " << GetRealTime() - real_normalization_start_time << "s.\n";
  PoreModel pore_model;
  pore_model.Load(pore_model_file_path_);
  SequenceBatch reference;
  reference.InitializeLoading(reference_file_path_);
  uint32_t num_reference_sequences = reference.LoadAllSequences();
  SignalBatch reference_signal_batch;
  reference_signal_batch.ConvertSequencesToSignals(reference, pore_model, num_reference_sequences);
  real_normalization_start_time = GetRealTime();
  for (size_t reference_signal_index = 0; reference_signal_index < num_reference_sequences; ++reference_signal_index) {
    reference_signal_batch.NormalizeSignalAt(reference_signal_index);
  }
  std::cerr << "Normalize " << num_reference_sequences << " reference signals in " << GetRealTime() - real_normalization_start_time << "s.\n";
  double real_start_time = GetRealTime();
  for (size_t read_signal_index = 0; read_signal_index < num_loaded_read_signals; ++read_signal_index) {
    for (size_t reference_signal_index = 0; reference_signal_index < num_reference_sequences; ++reference_signal_index) {
      std::cerr << "Read name: " << read_signal_batch.GetSignalNameAt(read_signal_index) << ", reference name: " << reference.GetSequenceNameAt(reference_signal_index) << "\n";
      sDTW(reference_signal_batch.GetSignalAt(reference_signal_index), read_signal_batch.GetSignalAt(read_signal_index));
    }
    std::cerr << num_loaded_read_signals << "\n";
  }
  std::cerr << "Finished mapping in " << GetRealTime() - real_start_time << ", # reads: " << num_loaded_read_signals << "\n";
  read_signal_batch.FinalizeLoading();
  reference.FinalizeLoading();
  reference_signal_batch.FinalizeLoading();
}

void Sigmap::CWTAlign() {
  SignalBatch read_signal_batch;
  read_signal_batch.InitializeLoading(signal_directory_);
  size_t num_loaded_read_signals = read_signal_batch.LoadAllReadSignals();
  PoreModel pore_model;
  pore_model.Load(pore_model_file_path_);
  SequenceBatch reference_sequence_batch;
  reference_sequence_batch.InitializeLoading(reference_file_path_);
  uint32_t num_reference_sequences = reference_sequence_batch.LoadAllSequences();
  SignalBatch reference_signal_batch;
  reference_signal_batch.ConvertSequencesToSignals(reference_sequence_batch, pore_model, num_reference_sequences);
  std::vector<std::vector<float> > reference_feature_signals;
  std::vector<std::vector<size_t> > reference_feature_positions;
  float cwt_scale0 = 1;
  for (size_t reference_signal_index = 0; reference_signal_index < num_reference_sequences; ++reference_signal_index) {
    reference_feature_signals.push_back(std::vector<float>());
    reference_feature_positions.push_back(std::vector<size_t>());
    GenerateFeatureSignalUsingCWT(reference_signal_batch.GetSignalAt(reference_signal_index), cwt_scale0, reference_feature_signals.back(), reference_feature_positions.back());
  }
  double real_start_time = GetRealTime();
  std::vector<float> read_feature_signal;
  std::vector<size_t> read_feature_positions;
  for (size_t read_signal_index = 0; read_signal_index < num_loaded_read_signals; ++read_signal_index) {
    read_feature_signal.clear();
    read_feature_positions.clear();
    ssize_t feature_mapping_end_position = -1;
    GenerateFeatureSignalUsingCWT(read_signal_batch.GetSignalAt(read_signal_index), 8 * cwt_scale0, read_feature_signal, read_feature_positions);
    for (size_t reference_signal_index = 0; reference_signal_index < num_reference_sequences; ++reference_signal_index) {
      std::cerr << "Read name: " << read_signal_batch.GetSignalNameAt(read_signal_index) << ", reference name: " << reference_sequence_batch.GetSequenceNameAt(reference_signal_index) << "\n";
      float dtw_distance = sDTW(reference_feature_signals[reference_signal_index].data(), reference_feature_signals[reference_signal_index].size(), read_feature_signal.data(), read_feature_signal.size(), feature_mapping_end_position);
      std::cerr << "DTW distance: " << dtw_distance << ", feature_mapping_end_position: " << feature_mapping_end_position << ", rough mapping end postion: " << reference_feature_positions[reference_signal_index][feature_mapping_end_position] << ".\n";
    }
    std::cerr << "\n";
  }
  std::cerr << "Finished mapping in " << GetRealTime() - real_start_time << "s, # reads: " << num_loaded_read_signals << "\n";
  read_signal_batch.FinalizeLoading();
  reference_sequence_batch.FinalizeLoading();
  reference_signal_batch.FinalizeLoading();
}

void Sigmap::ConstructIndex() {
  PoreModel pore_model;
  pore_model.Load(pore_model_file_path_);
  SequenceBatch reference_sequence_batch;
  reference_sequence_batch.InitializeLoading(reference_file_path_);
  uint32_t num_reference_sequences = reference_sequence_batch.LoadAllSequences();
  for (size_t reference_sequence_index = 0; reference_sequence_index < num_reference_sequences; ++reference_sequence_index) {
    reference_sequence_batch.PrepareNegativeSequenceAt(reference_sequence_index);
  }
  std::vector<std::vector<bool> > positive_is_masked;
  std::vector<std::vector<bool> > negative_is_masked;
  GenerateMaskedPositions(dimension_ + pore_model.GetKmerSize() - 1, 0.0002, num_reference_sequences, reference_sequence_batch, positive_is_masked, negative_is_masked);
  //
  //yak_ch_t *h;
  //yak_copt_t opt;
  //yak_copt_init(&opt);
  //opt.k = dimension_ + pore_model.GetKmerSize() - 1;
  //opt.n_thread = num_threads_;
  //h = yak_count_file(reference_file_path_.data(), NULL, &opt);
  //int64_t cnt[YAK_N_COUNTS];
  //yak_ch_hist(h, cnt, opt.n_thread);
  //yak_ch_get(h,);
  //
  SignalBatch reference_signal_batch;
  reference_signal_batch.ConvertSequencesToSignals(reference_sequence_batch, pore_model, num_reference_sequences);
  std::vector<std::vector<float> > positive_reference_feature_signals;
  std::vector<std::vector<float> > negative_reference_feature_signals;
  for (size_t reference_signal_index = 0; reference_signal_index < num_reference_sequences; ++reference_signal_index) {
    positive_reference_feature_signals.push_back(std::vector<float>());
    negative_reference_feature_signals.push_back(std::vector<float>());
    GenerateZscoreNormalizedSignal(reference_signal_batch.GetSignalAt(reference_signal_index).signal_values, reference_signal_batch.GetSignalLengthAt(reference_signal_index), positive_reference_feature_signals.back());
    GenerateZscoreNormalizedSignal(reference_signal_batch.GetSignalAt(reference_signal_index).negative_signal_values, reference_signal_batch.GetSignalLengthAt(reference_signal_index), negative_reference_feature_signals.back());
  }
  SpatialIndex spatial_index(dimension_, max_leaf_, 1, output_file_path_);
  spatial_index.Construct(positive_reference_feature_signals.size(), positive_is_masked, negative_is_masked, positive_reference_feature_signals, negative_reference_feature_signals);
  spatial_index.Save();
  reference_signal_batch.FinalizeLoading();
  //yak_ch_destroy(h);
}

void Sigmap::GenerateEvents(const Signal &signal, std::vector<float> &feature_signal) {
  std::vector<float> buffer;
  std::vector<Event> events;
  std::vector<float> prefix_sum;
  std::vector<float> prefix_sum_square;
  std::vector<float> tstat1;
  std::vector<float> tstat2;
  std::vector<size_t> peaks;
  const DetectorArgs ed_params = event_detection_defaults;
  DetectEvents(signal.signal_values, signal.signal_length, ed_params, prefix_sum, prefix_sum_square, tstat1, tstat2, peaks, events);
  feature_signal.clear();
  for (size_t ei = 0; ei < events.size(); ++ei) {
    feature_signal.emplace_back(events[ei].mean);
  }
  //feature_signal.swap(buffer);
  GenerateZscoreNormalizedSignal(feature_signal.data(), feature_signal.size(), buffer);
  feature_signal.clear();
  for (size_t i = 0; i < buffer.size(); i += 1) {
    //feature_signal.emplace_back(buffer[i]);
    //if (i >= 1) {
      if (i == 0 || (i >= 1 && (abs(buffer[i] - feature_signal.back()) > 0.1))) {
        feature_signal.emplace_back(buffer[i]);
      }
    //}
  }
#ifdef DEBUG
  std::cerr << "After compression: " << feature_signal.size() << "\n";
#endif
}

void Sigmap::GenerateFeatureSignalUsingCWT(const Signal &signal, float scale0, std::vector<float> &feature_signal, std::vector<size_t> &feature_positions) {
  std::vector<float> buffer;
  GenerateMADNormalizedSignal(signal.signal_values, signal.signal_length, buffer);
  GenerateCWTSignal(buffer.data(), buffer.size(), scale0, feature_signal);
  buffer.clear();
  float mean = GenerateZscoreNormalizedSignal(feature_signal.data(), feature_signal.size(), buffer);
  feature_signal.clear();
  GeneratePeaks(buffer.data(), buffer.size(), mean / 4, feature_signal, feature_positions);
}

float Sigmap::GenerateMADNormalizedSignal(const float *signal_values, size_t signal_length, std::vector<float> &normalized_signal) {
  // Should use a linear algorithm like median of medians
  // One such better algorithm can be found here: https://rcoh.me/posts/linear-time-median-finding/
  // But for now let us use sort
  normalized_signal.assign(signal_values, signal_values + signal_length);
  std::nth_element(normalized_signal.begin(), normalized_signal.begin() + signal_length / 2, normalized_signal.end());
  float signal_median = normalized_signal[signal_length / 2]; // This is a fake median, but should be okay for a quick implementation
  for (size_t i = 0; i < signal_length; ++i) {
    normalized_signal[i] = std::abs(normalized_signal[i] - signal_median);
  }
  std::nth_element(normalized_signal.begin(), normalized_signal.begin() + signal_length / 2, normalized_signal.end());
  float MAD = normalized_signal[signal_length / 2]; // Again, fake MAD, ok for a quick implementation
  // Now we can normalize signal
  for (size_t i = 0; i < signal_length; ++i) {
    normalized_signal[i] = (signal_values[i] - signal_median) / MAD;
  }
  return MAD;
}

float Sigmap::GenerateZscoreNormalizedSignal(const float *signal_values, size_t signal_length, std::vector<float> &normalized_signal) {
  // Calculate mean
  float mean = 0;
  for (size_t i = 0; i < signal_length; ++i) {
    mean += signal_values[i];
  }
  mean /= signal_length;
  // Calculate standard deviation
  float SD = 0;
  for (size_t i = 0; i < signal_length; ++i) {
    SD += (signal_values[i] - mean) * (signal_values[i] - mean);
  }
  SD /= (signal_length - 1);
  SD = sqrt(SD); 
  // Now we can normalize signal
  for (size_t i = 0; i < signal_length; ++i) {
    normalized_signal.emplace_back((signal_values[i] - mean) / SD);
  }
  return SD;
}

void Sigmap::GenerateCWTSignal(const float *signal_values, size_t signal_length, float scale0, std::vector<float> &cwt_signal) {
  char wave[] = "dog";
  char type[] = "pow";
  double param = 2.0;
  double dt = 1;
  double dj = 1; // Separation bewteen scales.
  int J = 1;
  int N = signal_length;
  cwt_object wt = cwt_init(wave, param, N, dt, J);
  setCWTScales(wt, scale0, dj, type, 2.0);
  cwt(wt, signal_values);
  for (size_t i = 0; i < signal_length; ++i) {
    cwt_signal.push_back(wt->output[i].re); 
  }
  //cwt_summary(wt);
  cwt_free(wt);
}

void Sigmap::GeneratePeaks(const float *signal_values, size_t signal_length, float selective, std::vector<float> &peaks, std::vector<size_t> &peak_positions) {
  float previous_valley = signal_values[0];
  float previous_peak = signal_values[0];
  for (size_t i = 1; i < signal_length - 1; ++i) {
    if (signal_values[i] > signal_values[i - 1] && signal_values[i] >= signal_values[i + 1] && signal_values[i] >= previous_valley + selective) {
      peaks.push_back(signal_values[i]);
      peak_positions.push_back(i);
      previous_peak = signal_values[i];
    } else if (signal_values[i] < signal_values[i - 1] && signal_values[i] <= signal_values[i + 1] && signal_values[i] <= previous_peak - selective) {
      peaks.push_back(signal_values[i]);
      peak_positions.push_back(i);
      previous_valley = signal_values[i];
    }
  }
}

void Sigmap::EventsToText() {
  SignalBatch read_signal_batch;
  read_signal_batch.InitializeLoading(signal_directory_);
  size_t num_loaded_read_signals = read_signal_batch.LoadAllReadSignals();
  std::vector<Event> events;
  const DetectorArgs ed_params = event_detection_defaults;
  FILE *output_file = fopen((output_file_path_ + "_event").c_str(), "w");
  assert(output_file != NULL);
  for (size_t i = 0; i < num_loaded_read_signals; ++i) {
    std::vector<float> prefix_sum;
    std::vector<float> prefix_sum_square;
    std::vector<float> tstat1;
    std::vector<float> tstat2;
    std::vector<size_t> peaks;
    const Signal &read_signal = read_signal_batch.GetSignalAt(i);
    events.clear();
    DetectEvents(read_signal.signal_values, read_signal.signal_length, ed_params, prefix_sum, prefix_sum_square, tstat1, tstat2, peaks, events);
    std::vector<float> buffer;
    std::vector<float> normalized_events;
    for (size_t ei = 0; ei < events.size(); ++ei) {
      buffer.emplace_back(events[ei].mean);
    }
    //GenerateMADNormalizedSignal(buffer.data(), buffer.size(), normalized_events);
    GenerateZscoreNormalizedSignal(buffer.data(), buffer.size(), normalized_events);
    //fprintf(output_file, "%s\t", read_signal.name);
    for (size_t ei = 0; ei < events.size(); ++ei) {
      fprintf(output_file, "%f\n", normalized_events[ei]);
    }
    //fprintf(output_file, "%f\n", normalized_events[events.size() - 1]);
  }
  fclose(output_file);
  read_signal_batch.FinalizeLoading();
}

void Sigmap::FAST5ToText() {
  SignalBatch read_signal_batch;
  read_signal_batch.InitializeLoading(signal_directory_);
  size_t num_loaded_read_signals = read_signal_batch.LoadAllReadSignals();
  FILE *output_file = fopen((output_file_path_ + "_fast5").c_str(), "w");
  assert(output_file != NULL);
  for (size_t i = 0; i < num_loaded_read_signals; ++i) {
    const Signal &read_signal = read_signal_batch.GetSignalAt(i);
    //fprintf(output_file, "%s\t", read_signal.name);
    for (size_t signal_position = 0; signal_position < read_signal.signal_length - 1; ++signal_position) {
      //fprintf(output_file, "%f\t", read_signal.signal_values[signal_position]);
      fprintf(output_file, "%f\n", read_signal.signal_values[signal_position]);
    }
    fprintf(output_file, "%f\n", read_signal.signal_values[read_signal.signal_length - 1]);
  }
  fclose(output_file);
  read_signal_batch.FinalizeLoading();
}

float Sigmap::sDTW(const float *target_signal_values, size_t target_length, const float *query_signal_values, size_t query_length, ssize_t &mapping_end_position) {
  double real_start_time = GetRealTime();
  float min_dtw_distance = std::numeric_limits<float>::max();
  mapping_end_position = -1;
  std::vector<float> previous_row(query_length + 1, std::numeric_limits<float>::max());
  previous_row[0] = 0;
  std::vector<float> current_row(query_length + 1);
  for (size_t target_position = 1; target_position <= target_length; ++target_position) {
    current_row[0] = 0;
    for (size_t query_position = 1; query_position <= query_length; ++query_position) {
      float cost = std::abs(target_signal_values[target_position - 1] - query_signal_values[query_position - 1]);
      current_row[query_position] = cost + std::min({previous_row[query_position - 1], previous_row[query_position], current_row[query_position - 1]});
    }
    if (current_row[query_length] < min_dtw_distance) {
      min_dtw_distance = current_row[query_length];
      mapping_end_position = target_position;
    }
    current_row.swap(previous_row);
  }
  std::cerr << "Finished sDTW in " << GetRealTime() - real_start_time << ", target length: " << target_length << ", query length: " << query_length << "\n";
  return min_dtw_distance;
}

float Sigmap::sDTW(const Signal &target_signal, const Signal &query_signal) {
  double real_start_time = GetRealTime();
  size_t query_length = query_signal.signal_length;
  size_t target_length = target_signal.signal_length;
  float min_cost = std::numeric_limits<float>::max();
  size_t mapping_end_position = 0;
  std::vector<float> previous_row(query_length + 1, std::numeric_limits<float>::max());
  previous_row[0] = 0;
  std::vector<float> current_row(query_length + 1);
  for (size_t target_position = 1; target_position <= target_length; ++target_position) {
    current_row[0] = 0;
    for (size_t query_position = 1; query_position <= query_length; ++query_position) {
      float cost = std::abs(target_signal.signal_values[target_position - 1] - query_signal.signal_values[query_position - 1]);
      current_row[query_position] = cost + std::min({previous_row[query_position - 1], previous_row[query_position], current_row[query_position - 1]});
    }
    if (current_row[query_length] < min_cost) {
      min_cost = current_row[query_length];
      mapping_end_position = target_position;
    }
    current_row.swap(previous_row);
  }
  std::cerr << "Finished sDTW in " << GetRealTime() - real_start_time << ", target length: " << target_length << ", query length: " << query_length << "\n";
  std::cerr << "DTW distance: " << min_cost << ", mapping_end_position: " << mapping_end_position << ".\n";
  return min_cost;
}

void SigmapDriver::ParseArgsAndRun(int argc, char *argv[]) {
  cxxopts::Options options("sigmap", "Map ONT raw signal data");
  options.add_options("Indexing")
    ("i,build-index", "Build spatial index for reference")
    ("d,dimension", "Dimension of spatial index", cxxopts::value<int>(), "INT")
    ("l,max-leaf", "Max leaf of spatial index", cxxopts::value<int>(), "INT");
  //options.add_options("Signal data indexing")
  //  ("build-sig-index", "Build index for signal data directory");
  options.add_options("Mapping")
    ("m,map", "Map signal data")
    ("t,num-threads", "# threads for mapping [1]", cxxopts::value<int>(), "INT");
  options.add_options("Input")
    ("r,ref", "Reference file", cxxopts::value<std::string>(), "FILE")
    ("p,pore-model", "Pore model file", cxxopts::value<std::string>(), "FILE")
    ("x,ref-index", "Reference index file", cxxopts::value<std::string>(), "FILE")
    //("sig-index", "Signal data directory index file", cxxopts::value<std::string>(), "FILE")
    ("s,sig-dir", "Signal data directory", cxxopts::value<std::string>(), "DIR");
    //("b,read-file", "Basecalled FASTA/FASTQ read file", cxxopts::value<std::string>());
  options.add_options("Output")
    ("o,output", "Output file", cxxopts::value<std::string>());
  options.add_options()
    ("h,help", "Print help");

  auto result = options.parse(argc, argv);
  int num_threads = 1;
  if (result.count("t")) {
    num_threads = result["num-threads"].as<int>();
  }

  if (result.count("i")) {
    int dimension = 6;
    if (result.count("d")) {
      dimension = result["dimension"].as<int>();
    }
    int max_leaf = 10;
    if (result.count("l")) {
      max_leaf = result["max-leaf"].as<int>();
    }
    std::cerr << "Dimension: " << dimension << ", max leaf: " << max_leaf << "\n";
    std::string reference_file_path;
    if (result.count("r")) {
      reference_file_path = result["ref"].as<std::string>();
    } else {
      sigmap::ExitWithMessage("No reference file specified!");
    }
    std::cerr << "Reference file: " << reference_file_path << "\n";
    std::string pore_model_file_path;
    if (result.count("p")) {
      pore_model_file_path = result["pore-model"].as<std::string>();
    } else {
      sigmap::ExitWithMessage("No pore model file specified!");
    }
    std::cerr << "Pore model file: " << pore_model_file_path << "\n";
    std::string output_file_path;
    if (result.count("o")) {
      output_file_path = result["output"].as<std::string>();
    } else {
      sigmap::ExitWithMessage("No output file specified!");
    }
    std::cerr << "Output file: " << output_file_path << "\n";
    Sigmap sigmap_for_indexing(dimension, max_leaf, reference_file_path, pore_model_file_path, output_file_path);
    sigmap_for_indexing.ConstructIndex();
  } else if (result.count("m")) {
    std::cerr << "Number of threads: " << num_threads << "\n";
    std::string reference_file_path;
    if (result.count("r")) {
      reference_file_path = result["ref"].as<std::string>();
    } else {
      sigmap::ExitWithMessage("No reference file specified!");
    }
    std::cerr << "Reference file: " << reference_file_path << "\n";
    std::string pore_model_file_path;
    if (result.count("p")) {
      pore_model_file_path = result["pore-model"].as<std::string>();
    } else {
      sigmap::ExitWithMessage("No pore model file specified!");
    }
    std::cerr << "Pore model file: " << pore_model_file_path << "\n";
    std::string reference_index_file_path;
    if (result.count("x")) {
      reference_index_file_path = result["ref-index"].as<std::string>();
    } else {
      sigmap::ExitWithMessage("No reference index file specified!");
    }
    std::cerr << "Reference index file: " << reference_index_file_path << "\n";
    std::string signal_dir;
    if (result.count("sig-dir")) {
      signal_dir = result["sig-dir"].as<std::string>();
    } else {
      sigmap::ExitWithMessage("No signal data directory specified!");
    }
    std::cerr << "Signal directory: " << signal_dir << "\n";
    std::string output_file_path;
    if (result.count("o")) {
      output_file_path = result["output"].as<std::string>();
    } else {
      sigmap::ExitWithMessage("No output file specified!");
    }
    std::cerr << "Output file: " << output_file_path << "\n";
    Sigmap sigmap_for_mapping(num_threads, reference_file_path, pore_model_file_path, signal_dir, reference_index_file_path, output_file_path);
    //sigmap_for_mapping.CWTAlign();
    //sigmap_for_mapping.DTWAlign();
    sigmap_for_mapping.Map();
    //sigmap_for_mapping.FAST5ToText();
    //sigmap_for_mapping.EventsToText();
  } else if (result.count("h")) {
    std::cerr << options.help({"", "Indexing", "Mapping", "Input", "Output"});
  } else {
    std::cerr << options.help({"", "Indexing", "Mapping", "Input", "Output"});
  }
  //std::string read_file_path;
  //if (result.count("b")) {
  //  read_file_path = result["read-file"].as<std::string>();
  //} else {
  //  sigmap::ExitWithMessage("No read file specified!");
  //}
  //std::cerr << "Read file: " << read_file_path << "\n";
}
} // namespace sigmap

int main(int argc, char *argv[]) {
  sigmap::SigmapDriver sigmap_driver;
  sigmap_driver.ParseArgsAndRun(argc, argv);
  return 0;
}
