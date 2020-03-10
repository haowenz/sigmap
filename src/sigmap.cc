#include "sigmap.h"

#include <cassert>
#include <iostream>
#include <string>

#include "cxxopts.hpp"
#include "sequence_batch.h"
#include "pore_model.h"
#include "utils.h"

namespace sigmap {
void Sigmap::Map() {
  SignalBatch read_signal_batch;
  read_signal_batch.InitializeLoading(signal_directory_);
  size_t num_loaded_read_signals = read_signal_batch.LoadAllReadSignals();
  //double real_normalization_start_time = GetRealTime();
  //for (size_t read_index = 0; read_index < num_loaded_read_signals; ++read_index) {
  //  read_signal_batch.NormalizeSignalAt(read_index);
  //}
  //std::cerr << "Normalize " << num_loaded_read_signals << " read signals in " << GetRealTime() - real_normalization_start_time << "s.\n";
  PoreModel pore_model;
  pore_model.Load(pore_model_file_path_);
  SequenceBatch reference;
  reference.InitializeLoading(reference_file_path_);
  uint32_t num_reference_sequences = reference.LoadAllSequences();
  SignalBatch reference_signal_batch;
  reference_signal_batch.ConvertSequencesToSignals(reference, pore_model, num_reference_sequences);
  //real_normalization_start_time = GetRealTime();
  //for (size_t reference_signal_index = 0; reference_signal_index < num_reference_sequences; ++reference_signal_index) {
  //  reference_signal_batch.NormalizeSignalAt(reference_signal_index);
  //}
  //std::cerr << "Normalize " << num_reference_sequences << " reference signals in " << GetRealTime() - real_normalization_start_time << "s.\n";
  //for (size_t read_signal_index = 0; read_signal_index < num_loaded_read_signals; ++read_signal_index) {
  //  for (size_t reference_signal_index = 0; reference_signal_index < num_reference_sequences; ++reference_signal_index) {
  //    std::cerr << "Read name: " << read_signal_batch.GetSignalNameAt(read_signal_index) << ", reference name: " << reference.GetSequenceNameAt(reference_signal_index) << "\n";
  //    sDTW(reference_signal_batch.GetSignalAt(reference_signal_index), read_signal_batch.GetSignalAt(read_signal_index));
  //  }
  //}
  std::vector<std::vector<float> > reference_signals;
  for (size_t reference_signal_index = 0; reference_signal_index < num_reference_sequences; ++reference_signal_index) {
    reference_signals.push_back(std::vector<float>());
    GetFeatureSignal(reference_signal_batch.GetSignalAt(reference_signal_index), sqrt(2), reference_signals.back());
  }
  std::vector<float> read_signal;
  for (size_t read_signal_index = 0; read_signal_index < num_loaded_read_signals; ++read_signal_index) {
    read_signal.clear();
    GetFeatureSignal(read_signal_batch.GetSignalAt(read_signal_index), 8 * sqrt(2), read_signal);
    for (size_t reference_signal_index = 0; reference_signal_index < num_reference_sequences; ++reference_signal_index) {
      std::cerr << "Read name: " << read_signal_batch.GetSignalNameAt(read_signal_index) << ", reference name: " << reference.GetSequenceNameAt(reference_signal_index) << "\n";
      sDTW(reference_signals[reference_signal_index].data(), reference_signals[reference_signal_index].size(), read_signal.data(), read_signal.size());
    }
    std::cerr << "\n";
  }
  read_signal_batch.FinalizeLoading();
  reference.FinalizeLoading();
  reference_signal_batch.FinalizeLoading();
}

void Sigmap::GetFeatureSignal(const Signal &signal, float scale0, std::vector<float> &signal_feature) {
  std::vector<float> buffer;
  GetNormalizedSignal(signal.signal, signal.signal_length, buffer);
  GetCWTSignal(buffer.data(), buffer.size(), scale0, signal_feature);
  buffer.clear();
  float MAD = GetNormalizedSignal(signal_feature.data(), signal_feature.size(), buffer);
  signal_feature.clear();
  GetPeaks(buffer.data(), buffer.size(), MAD / 8, signal_feature);
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

float Sigmap::sDTW(const float *target_signal, size_t target_length, const float *query_signal, size_t query_length) {
  double real_start_time = GetRealTime();
  //size_t query_length = query_signal.signal_length;
  //size_t target_length = target_signal.signal_length;
  float min_dtw_distance = std::numeric_limits<float>::max();
  size_t mapping_end_position = 0;
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
  std::cerr << "DTW distance: " << min_dtw_distance << ", mapping_end_position: " << mapping_end_position << ".\n";
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
    ("i,build-index", "Build reference index");
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
    Sigmap sigmap_for_mapping(reference_file_path, pore_model_file_path, signal_dir, output_file_path);
    sigmap_for_mapping.Map();
    //sigmap_for_mapping.FAST5ToText();
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
