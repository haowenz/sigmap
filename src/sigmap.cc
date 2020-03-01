#include "sigmap.h"

#include <cassert>
#include <iostream>
#include <string>

#include "cxxopts.hpp"
#include "pore_model.h"
#include "signal_batch.h"
#include "utils.h"

namespace sigmap {
void Sigmap::Map() {
  SignalBatch signal_batch;
  signal_batch.InitializeLoading(signal_directory_);
  size_t num_loaded_read_signals = signal_batch.LoadAllReadSignals();
  signal_batch.FinalizeLoading();
  PoreModel pore_model;
  pore_model.Load(pore_model_file_path_);
  double real_normalization_start_time = GetRealTime();
  for (size_t read_index = 0; read_index < num_loaded_read_signals; ++read_index) {
    signal_batch.NormalizeSignalAt(read_index);
  }
  std::cerr << "Normalize " << num_loaded_read_signals << " read signals in " << GetRealTime() - real_normalization_start_time << "s.\n";
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
      reference_file_path = result["ref-index"].as<std::string>();
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
    Sigmap sigmap_for_mapping(pore_model_file_path, signal_dir);
    sigmap_for_mapping.Map();
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
