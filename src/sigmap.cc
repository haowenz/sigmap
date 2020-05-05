#include "sigmap.h"

#include <cassert>
#include <iostream>
#include <string>

#include "cxxopts.hpp"
#include "nanoflann.hpp"
#include "spatial_index.h"
#include "sequence_batch.h"
#include "pore_model.h"
#include "utils.h"

namespace sigmap {
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
  // Use pore model to convert reference sequence to signal
  SignalBatch reference_signal_batch;
  reference_signal_batch.ConvertSequencesToSignals(reference_sequence_batch, pore_model, num_reference_sequences);
  // Perform CWT on reference signals
  std::vector<std::vector<float> > reference_feature_signals;
  std::vector<std::vector<size_t> > reference_feature_positions;
  float cwt_scale0 = 1;
  for (size_t reference_signal_index = 0; reference_signal_index < num_reference_sequences; ++reference_signal_index) {
    reference_feature_signals.push_back(std::vector<float>());
    reference_feature_positions.push_back(std::vector<size_t>());
    GetFeatureSignal(reference_signal_batch.GetSignalAt(reference_signal_index), cwt_scale0, reference_feature_signals.back(), reference_feature_positions.back());
  }
  // Load spatial index for reference signals 
  SpatialIndex reference_spatial_index(1000, std::vector<int>(1000,5000), reference_index_file_path_);
  reference_spatial_index.Load();
  // Map each reads
  double real_start_time = GetRealTime();
  int read_signal_point_cloud_step_size = 1;
  float search_radius = 0.3;
  std::vector<float> read_feature_signal;
  std::vector<size_t> read_feature_positions;
  std::vector<std::vector<float> > read_point_cloud;
  std::vector<uint64_t> positive_hits;
  std::vector<uint64_t> negative_hits;
  std::vector<uint64_t> positive_candidates;
  std::vector<uint64_t> negative_candidates;
  std::vector<SignalAnchorChain> positive_chains;
  for (size_t read_signal_index = 0; read_signal_index < num_loaded_read_signals; ++read_signal_index) {
    read_feature_signal.clear();
    read_feature_positions.clear();
    ssize_t feature_mapping_end_position = -1;
    GetReadFeatureSignal(read_signal_batch.GetSignalAt(read_signal_index), cwt_scale0, read_feature_signal, read_feature_positions);
    read_point_cloud.clear();
    positive_hits.clear();
    negative_hits.clear();
    positive_candidates.clear();
    negative_candidates.clear();
    positive_chains.clear();
    reference_spatial_index.GeneratePointCloud(read_feature_signal.data(), read_feature_signal.size(), read_signal_point_cloud_step_size, read_point_cloud);
    //reference_spatial_index.GenerateCandidates(read_point_cloud, &positive_hits, &negative_hits, &positive_candidates, &negative_candidates);
    reference_spatial_index.GenerateChains(read_point_cloud, read_signal_point_cloud_step_size, search_radius, num_reference_sequences, reference_feature_signals, positive_chains);
    std::cerr << "Max chaining score: " << positive_chains[0].score << ", signal_index: " << positive_chains[0].start_position << ", anchor target start postion: " << positive_chains[0].end_position << ", rough end position on ref: " << reference_feature_positions[positive_chains[0].start_position][positive_chains[0].end_position] << ".\n";
    std::cerr << "Read name: " << read_signal_batch.GetSignalNameAt(read_signal_index) << ", length: " << read_feature_signal.size() << ", reference name: " << reference_sequence_batch.GetSequenceNameAt(positive_chains[0].start_position) << ", length: " << reference_feature_signals[positive_chains[0].start_position].size() << "\n";
    std::cerr << "\n";
  }
  std::cerr << "Finished mapping in " << GetRealTime() - real_start_time << ", # reads: " << num_loaded_read_signals << "\n";
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
    GetFeatureSignal(reference_signal_batch.GetSignalAt(reference_signal_index), cwt_scale0, reference_feature_signals.back(), reference_feature_positions.back());
  }
  double real_start_time = GetRealTime();
  std::vector<float> read_feature_signal;
  std::vector<size_t> read_feature_positions;
  for (size_t read_signal_index = 0; read_signal_index < num_loaded_read_signals; ++read_signal_index) {
    read_feature_signal.clear();
    read_feature_positions.clear();
    ssize_t feature_mapping_end_position = -1;
    GetReadFeatureSignal(read_signal_batch.GetSignalAt(read_signal_index), 8 * cwt_scale0, read_feature_signal, read_feature_positions);
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
  SignalBatch reference_signal_batch;
  reference_signal_batch.ConvertSequencesToSignals(reference_sequence_batch, pore_model, num_reference_sequences);
  std::vector<std::vector<float> > reference_feature_signals;
  std::vector<std::vector<size_t> > reference_feature_positions;
  float cwt_scale0 = 1;
  for (size_t reference_signal_index = 0; reference_signal_index < num_reference_sequences; ++reference_signal_index) {
    reference_feature_signals.push_back(std::vector<float>());
    reference_feature_positions.push_back(std::vector<size_t>());
    GetFeatureSignal(reference_signal_batch.GetSignalAt(reference_signal_index), cwt_scale0, reference_feature_signals.back(), reference_feature_positions.back());
  }
  SpatialIndex spatial_index(dimension_, max_leaf_, 1, output_file_path_);
  spatial_index.Construct(reference_feature_signals.size(), reference_feature_signals);
  spatial_index.Save();
  reference_signal_batch.FinalizeLoading();
}

void Sigmap::GetFeatureSignal(const Signal &signal, float scale0, std::vector<float> &feature_signal, std::vector<size_t> &feature_positions) {
  ////SaveVectorToFile(signal.signal, signal.signal_length, output_file_path_ + "_" + std::string(signal.name) + ".raw");
  //std::vector<float> buffer;
  //GetNormalizedSignal(signal.signal, signal.signal_length, buffer);
  ////SaveVectorToFile(buffer.data(), buffer.size(), output_file_path_ + "_" + std::string(signal.name) + ".normalized");
  //GetCWTSignal(buffer.data(), buffer.size(), scale0, feature_signal);
  ////SaveVectorToFile(feature_signal.data(), feature_signal.size(), output_file_path_ + "_" + std::string(signal.name) + ".cwt");
  //buffer.clear();
  ////float MAD = GetNormalizedSignal(feature_signal.data(), feature_signal.size(), buffer);
  //float MAD = GetZscoreNormalizedSignal(feature_signal.data(), feature_signal.size(), buffer);
  ////SaveVectorToFile(buffer.data(), buffer.size(), output_file_path_ + "_" + std::string(signal.name) + ".normalized_cwt");
  ////buffer.swap(feature_signal);
  //feature_signal.clear();
  //GetPeaks(buffer.data(), buffer.size(), MAD / 4, feature_signal, feature_positions);
  ////SaveVectorToFile(feature_signal.data(), feature_signal.size(), output_file_path_ + "_" + std::string(signal.name) + ".peaks");
  ////SaveVectorToFile(feature_positions.data(), feature_positions.size(), output_file_path_ + "_" + std::string(signal.name) + ".peak_positions");
  GetNormalizedSignal(signal.signal, signal.signal_length, feature_signal);
  for (size_t i = 0; i < signal.signal_length; ++i) {
    feature_positions.emplace_back(i);
  }
}

void Sigmap::GetReadFeatureSignal(const Signal &signal, float scale0, std::vector<float> &feature_signal, std::vector<size_t> &feature_positions) {
  //std::vector<Event> events;
  //const DetectorArgs ed_params = event_detection_defaults;
  //detect_events(signal.signal, signal.signal_length, ed_params, events);
  //for (size_t ei = 0; ei < events.size(); ++ei) {
  //  feature_signal.emplace_back(events[ei].mean);
  //}
  ////SaveVectorToFile(signal.signal, signal.signal_length, output_file_path_ + "_" + std::string(signal.name) + ".raw");
  //std::vector<float> buffer;
  //GetNormalizedSignal(feature_signal.data(), feature_signal.size(), buffer);
  //feature_signal.clear();
  ////SaveVectorToFile(buffer.data(), buffer.size(), output_file_path_ + "_" + std::string(signal.name) + ".normalized");
  //GetCWTSignal(buffer.data(), buffer.size(), scale0, feature_signal);
  ////SaveVectorToFile(feature_signal.data(), feature_signal.size(), output_file_path_ + "_" + std::string(signal.name) + ".cwt");
  //buffer.clear();
  ////float MAD = GetNormalizedSignal(feature_signal.data(), feature_signal.size(), buffer);
  //float MAD = GetZscoreNormalizedSignal(feature_signal.data(), feature_signal.size(), buffer);
  ////SaveVectorToFile(buffer.data(), buffer.size(), output_file_path_ + "_" + std::string(signal.name) + ".normalized_cwt");
  ////buffer.swap(feature_signal);
  //feature_signal.clear();
  //GetPeaks(buffer.data(), buffer.size(), MAD / 4, feature_signal, feature_positions);
  ////SaveVectorToFile(feature_signal.data(), feature_signal.size(), output_file_path_ + "_" + std::string(signal.name) + ".peaks");
  ////SaveVectorToFile(feature_positions.data(), feature_positions.size(), output_file_path_ + "_" + std::string(signal.name) + ".peak_positions");

  std::vector<float> buffer;
  std::vector<Event> events;
  const DetectorArgs ed_params = event_detection_defaults;
  detect_events(signal.signal, signal.signal_length, ed_params, events);
  //for (size_t ei = 0; ei < events.size(); ++ei) {
  //  buffer.emplace_back(events[ei].mean);
  //}
  //GetNormalizedSignal(buffer.data(), buffer.size(), feature_signal);
  //for (size_t i = 0; i < buffer.size(); ++i) {
  //  feature_positions.emplace_back(i);
  //}
  feature_signal.clear();
  for (size_t ei = 0; ei < events.size(); ++ei) {
    feature_signal.emplace_back(events[ei].mean);
  }
  GetNormalizedSignal(feature_signal.data(), feature_signal.size(), buffer);
  feature_signal.clear();
  for (size_t i = 0; i < buffer.size(); i += 2) {
    feature_signal.emplace_back(buffer[i]);
    feature_positions.emplace_back(i);
  }

  //std::vector<float> buffer;
  //GetNormalizedSignal(signal.signal, signal.signal_length, buffer);
  //std::vector<Event> events;
  //const DetectorArgs ed_params = event_detection_defaults;
  //detect_events(buffer.data(), buffer.size(), ed_params, events);
  //feature_signal.clear();
  //for (size_t i = 0; i < buffer.size(); i += 2) {
  //  feature_signal.emplace_back(buffer[i]);
  //  feature_positions.emplace_back(i);
  //}
}

void Sigmap::EventsToText() {
  SignalBatch read_signal_batch;
  read_signal_batch.InitializeLoading(signal_directory_);
  size_t num_loaded_read_signals = read_signal_batch.LoadAllReadSignals();
  std::vector<Event> events;
  const DetectorArgs ed_params = event_detection_defaults;
  FILE *output_file = fopen(output_file_path_.c_str(), "w");
  assert(output_file != NULL);
  for (size_t i = 0; i < num_loaded_read_signals; ++i) {
    const Signal &read_signal = read_signal_batch.GetSignalAt(i);
    events.clear();
    detect_events(read_signal.signal, read_signal.signal_length, ed_params, events);
    std::vector<float> buffer;
    std::vector<float> normalized_events;
    for (size_t ei = 0; ei < events.size(); ++ei) {
      buffer.emplace_back(events[ei].mean);
    }
    GetNormalizedSignal(buffer.data(), buffer.size(), normalized_events);
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
  FILE *output_file = fopen(output_file_path_.c_str(), "w");
  assert(output_file != NULL);
  for (size_t i = 0; i < num_loaded_read_signals; ++i) {
    const Signal &read_signal = read_signal_batch.GetSignalAt(i);
    fprintf(output_file, "%s\t", read_signal.name);
    for (size_t signal_position = 0; signal_position < read_signal.signal_length - 1; ++signal_position) {
      fprintf(output_file, "%f\t", read_signal.signal[signal_position]);
    }
    fprintf(output_file, "%f\n", read_signal.signal[read_signal.signal_length - 1]);
  }
  fclose(output_file);
  read_signal_batch.FinalizeLoading();
}

float Sigmap::sDTW(const float *target_signal, size_t target_length, const float *query_signal, size_t query_length, ssize_t &mapping_end_position) {
  double real_start_time = GetRealTime();
  float min_dtw_distance = std::numeric_limits<float>::max();
  mapping_end_position = -1;
  std::vector<float> previous_row(query_length + 1, std::numeric_limits<float>::max());
  previous_row[0] = 0;
  std::vector<float> current_row(query_length + 1);
  for (size_t target_position = 1; target_position <= target_length; ++target_position) {
    current_row[0] = 0;
    for (size_t query_position = 1; query_position <= query_length; ++query_position) {
      float cost = std::abs(target_signal[target_position - 1] - query_signal[query_position - 1]);
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
      float cost = std::abs(target_signal.signal[target_position - 1] - query_signal.signal[query_position - 1]);
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
    int dimension = 5;
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
    Sigmap sigmap_for_mapping(reference_file_path, pore_model_file_path, signal_dir, reference_index_file_path, output_file_path);
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
