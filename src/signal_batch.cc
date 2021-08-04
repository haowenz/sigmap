#include "signal_batch.h"

#include <hdf5.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "utils.h"

namespace sigmap {
void SignalBatch::InitializeLoading(const std::string &signal_directory) {
  if (IsDirectory(signal_directory)) {
    signal_directory_ = signal_directory;
  } else {
    ExitWithMessage("Signal directory is in valid!");
  }
}

void SignalBatch::FinalizeLoading() {}

size_t SignalBatch::LoadAllReadSignals() {
  double real_start_time = GetRealTime();
  auto dir_list = ListDirectory(signal_directory_);
  // two rounds, get all the fast5 absolute paths first and then load reads
  std::vector<std::string> fast5_list;
  for (size_t pi = 0; pi < dir_list.size(); ++pi) {
    std::string &relative_path = dir_list[pi];
    if (relative_path == "." or relative_path == "..") {
      continue;
    }
    std::string absolute_path = signal_directory_ + "/" + relative_path;
    if (IsDirectory(absolute_path)) {
      auto sub_directory_list = ListDirectory(absolute_path);
      std::vector<std::string> sub_relative_paths;
      sub_relative_paths.reserve(sub_directory_list.size());
      for (size_t si = 0; si < sub_directory_list.size(); ++si) {
        sub_relative_paths.push_back(relative_path + "/" +
                                     sub_directory_list[si]);
      }
      dir_list.insert(dir_list.end(), sub_relative_paths.begin(),
                      sub_relative_paths.end());
    }
    bool is_fast5 = absolute_path.find(".fast5") != std::string::npos;
    if (is_fast5) {
      fast5_list.emplace_back(absolute_path);
    } else {
      // don't put it into index
    }
  }
  signals_.reserve(fast5_list.size());
  for (size_t pi = 0; pi < fast5_list.size(); ++pi) {
    AddSignalsFromFAST5(fast5_list[pi]);
  }
  std::cerr << "Loaded " << signals_.size() << " reads in "
            << GetRealTime() - real_start_time << "s.\n";
  return signals_.size();
}

void SignalBatch::AddSignalsFromFAST5(const std::string &fast5_file_path) {
  hdf5_tools::File fast5_file;
  fast5_file.open(fast5_file_path);
  // Check format
  bool is_single = false;
  std::vector<std::string> fast5_file_groups = fast5_file.list_group("/");
  for (std::string &group : fast5_file_groups) {
    if (group == "Raw") {
      is_single = true;
      break;
    }
  }
  if (is_single) {
    for (std::string &read : fast5_file.list_group("/Raw/Reads")) {
      std::string read_id = "";
      for (auto a : fast5_file.get_attr_map("/Raw/Reads/" + read)) {
        if (a.first == "read_id") {
          read_id = a.second;
          break;
        }
      }
      if (read_id.empty()) {
        std::cerr << "Error: failed to find read_id\n";
        exit(-1);
      }
      std::string raw_path = "/Raw/Reads/" + read;
      std::string ch_path = "/UniqueGlobalKey/channel_id";
      AddSignal(fast5_file, raw_path, ch_path);
    }
  } else {
    signals_.reserve(signals_.size() + fast5_file_groups.size());
    for (std::string &read : fast5_file_groups) {
      std::string raw_path = "/" + read + "/Raw";
      std::string ch_path = "/" + read + "/channel_id";
      AddSignal(fast5_file, raw_path, ch_path);
    }
  }
  fast5_file.close();
}

void SignalBatch::AddSignal(const hdf5_tools::File &file,
                            const std::string &raw_path,
                            const std::string &ch_path) {
  std::string id;
  for (auto a : file.get_attr_map(raw_path)) {
    if (a.first == "read_id") {
      id = a.second;
    } else if (a.first == "read_number") {
      // number = atoi(a.second.c_str());
    } else if (a.first == "start_time") {
      // start_sample = atoi(a.second.c_str());
    }
  }

  float digitisation = 0, range = 0, offset = 0;
  for (auto a : file.get_attr_map(ch_path)) {
    if (a.first == "channel_number") {
      // channel_idx = atoi(a.second.c_str()) - 1;
    } else if (a.first == "digitisation") {
      digitisation = atof(a.second.c_str());
    } else if (a.first == "range") {
      range = atof(a.second.c_str());
    } else if (a.first == "offset") {
      offset = atof(a.second.c_str());
    }
  }

  std::string sig_path = raw_path + "/Signal";
  std::vector<float> signal_values;
  file.read(sig_path, signal_values);
  // convert to pA
  uint32_t valid_signal_length = 0;
  float scale = range / digitisation;
  for (size_t i = 0; i < signal_values.size(); i++) {
    if ((signal_values[i] + offset) * scale > 30 &&
        (signal_values[i] + offset) * scale < 200) {
      signal_values[valid_signal_length] = (signal_values[i] + offset) * scale;
      ++valid_signal_length;
    }
  }
  if (valid_signal_length < signal_values.size()) {
    signal_values.erase(signal_values.begin() + valid_signal_length,
                        signal_values.end());
  }
  signals_.emplace_back(Signal{id, digitisation, range, offset, signal_values,
                               std::vector<float>()});
}

void SignalBatch::NormalizeSignalAt(size_t signal_index) {
  //// Should use a linear algorithm like median of medians
  //// But for now let us use sort
  // std::vector<float> tmp_signal((signals_[signal_index]).signal_values,
  // (signals_[signal_index]).signal_values +
  // (signals_[signal_index]).signal_length);
  // std::nth_element(tmp_signal.begin(), tmp_signal.begin() + tmp_signal.size()
  // / 2, tmp_signal.end()); float signal_median = tmp_signal[tmp_signal.size()
  // / 2]; // This is a fake median, but should be okay for a quick
  // implementation for (size_t i = 0; i < tmp_signal.size(); ++i) {
  //  tmp_signal[i] = std::abs(tmp_signal[i] - signal_median);
  //}
  // std::nth_element(tmp_signal.begin(), tmp_signal.begin() + tmp_signal.size()
  // / 2, tmp_signal.end()); float MAD = tmp_signal[tmp_signal.size() / 2]; //
  // Again, fake MAD, ok for a quick implementation
  //// Now we can normalize signal
  // for (size_t i = 0; i < signals_[signal_index].signal_length; ++i) {
  //  ((signals_[signal_index]).signal_values)[i] =
  //  (((signals_[signal_index]).signal_values)[i] - signal_median) / MAD;
  //}
  // Calculate mean
  float mean = 0;
  for (size_t i = 0; i < signals_[signal_index].GetSignalLength(); ++i) {
    mean += signals_[signal_index].signal_values[i];
  }
  mean /= signals_[signal_index].signal_values.size();
  // Calculate standard deviation
  float stdv = 0;
  for (size_t i = 0; i < signals_[signal_index].GetSignalLength(); ++i) {
    stdv += (signals_[signal_index].signal_values[i] - mean) *
            (signals_[signal_index].signal_values[i] - mean);
  }
  stdv /= (signals_[signal_index].signal_values.size() - 1);
  stdv = sqrt(stdv);
  // Now we can normalize signal
  for (size_t i = 0; i < signals_[signal_index].signal_values.size(); ++i) {
    signals_[signal_index].signal_values[i] =
        (signals_[signal_index].signal_values[i] - mean) / stdv;
  }
}

void SignalBatch::ConvertSequencesToSignals(const SequenceBatch &sequence_batch,
                                            const PoreModel &pore_model,
                                            size_t num_sequences) {
  double real_start_time = GetRealTime();
  for (size_t sequence_index = 0; sequence_index < num_sequences;
       ++sequence_index) {
    size_t sequence_length = sequence_batch.GetSequenceLengthAt(sequence_index);
    std::vector<float> signal_values = pore_model.GetLevelMeansAt(
        sequence_batch.GetSequenceAt(sequence_index), 0, sequence_length);
    std::vector<float> negative_signal_values = pore_model.GetLevelMeansAt(
        sequence_batch.GetNegativeSequenceAt(sequence_index).data(), 0,
        sequence_length);
    std::string signal_id(sequence_batch.GetSequenceNameAt(sequence_index));
    signals_.emplace_back(
        Signal{signal_id, 0, 0, 0, signal_values, negative_signal_values});
  }
  std::cerr << "Convert " << num_sequences << " sequences to signals in "
            << GetRealTime() - real_start_time << "s.\n";
}
}  // namespace sigmap
