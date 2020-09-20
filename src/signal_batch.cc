#include "signal_batch.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <hdf5.h>
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

void SignalBatch::FinalizeLoading() {
}

size_t SignalBatch::LoadAllReadSignals() {
  double real_start_time = GetRealTime();
  size_t num_reads = 0;
  auto dir_list = ListDirectory(signal_directory_);
  for (const auto& relative_path : dir_list) {
    if (relative_path == "." or relative_path == "..") {
      continue;
    }
    std::string absolute_path = signal_directory_ + "/" + relative_path;
    bool is_fast5 = absolute_path.find(".fast5") != std::string::npos;
    bool in_map = false;
    //bool in_map = fast5_to_read_name_map.find(fn) != fast5_to_read_name_map.end();
    // JTS 04/19: is_directory is painfully slow so we first check if the file is in the name map
    // if it is, it is definitely not a directory so we can skip the system call
    //if(!in_map && IsDirectory(absolute_path)) {
    if(in_map) {
      // recurse
      //index_path(read_db, full_fn, fast5_to_read_name_map);
    } else if (is_fast5) {
      if (in_map) {
        //index_file_from_map(read_db, full_fn, fast5_to_read_name_map);
      } else {
        AddSignalsFromFAST5(absolute_path);
        ++num_reads;
      }
    } else {
      // don't put it into index
    }
  }
  std::cerr << "Loaded " << num_reads << " reads in " << GetRealTime() - real_start_time << "s.\n";
  return num_reads;
}

void SignalBatch::AddSignalsFromFAST5(const std::string &fast5_file_path) {
  FAST5File fast5_file = OpenFAST5(fast5_file_path);
  if (fast5_file.is_multi_fast5) {
    ExitWithMessage("Multi-fast5 is not supported yet!");
    //std::vector<std::string> &read_names = GetReadNamesFromMultiFAST5(fast5_file);
    //std::string prefix = "read_";
    //for (size_t read_index = 0; read_index < read_names.size(); ++read_index) {
    //  std::string &read_name = read_names[read_index];
    //  if (read_name.find(prefix) == 0) {
    //    //std::string read_id = group_name.substr(prefix.size());
    //    //AddReadSignal(fast5_file, read_name);
    //  }
    //}
  } else {
    //std::string read_name = GetReadNameFromSingleFAST5(fast5_file);
    //AddReadSignal(fast5_file, read_name);
    AddSignalFromSingleFAST5(fast5_file);
  }
  CloseFAST5(fast5_file);
}

void SignalBatch::AddSignalFromSingleFAST5(const FAST5File& fast5_file) {
  //Signal read_signal = new Signal();
  size_t group_name_length;
  char *group_name = nullptr;
  //size_t read_name_length;
  char *read_name = nullptr;
  size_t signal_length;
  float *signal_values;
  float digitisation;
  float range;
  float offset;
  herr_t fast5_err;
  float scale;
  int num_dims;
  // Get read name length and name
  // Note that this weird function returns string length when setting buffer to null.
  // But when reading the name, the buffer size should be set to string length + 1!
  ssize_t read_group_name_length = H5Lget_name_by_idx(fast5_file.hdf5_file, LEGACY_FAST5_RAW_ROOT, H5_INDEX_NAME, H5_ITER_INC, 0, NULL, 0, H5P_DEFAULT);
  if (read_group_name_length < 0) {
    ExitWithMessage("The read name length is invalid!\n");
  }
  group_name_length = (size_t)read_group_name_length;
  group_name = (char*)calloc(1 + group_name_length, sizeof(char));
  read_group_name_length = H5Lget_name_by_idx(fast5_file.hdf5_file, LEGACY_FAST5_RAW_ROOT, H5_INDEX_NAME, H5_ITER_INC, 0, group_name, 1 + group_name_length, H5P_DEFAULT);
  if (read_group_name_length != (ssize_t)group_name_length) {
    ExitWithMessage("Read name lengths don't match! Failed to load read name.");
  }
  // Get read id and use it as read name
  std::string read_group = std::string(LEGACY_FAST5_RAW_ROOT) + "/" + std::string(group_name);
  hid_t read_group_id = H5Gopen(fast5_file.hdf5_file, read_group.data(), H5P_DEFAULT);
  if (read_group_id < 0) {
    fprintf(stderr, "Failed to open read group %s\n", read_group.data());
    exit(-1); // TODO(Haowen): fix this later
  }
  GetStringAttributeInGroup(read_group_id, "read_id", &read_name);
  // Get channel parameters
  const char *channel_id_group = "/UniqueGlobalKey/channel_id";
  hid_t channel_id_group_id = H5Gopen(fast5_file.hdf5_file, channel_id_group, H5P_DEFAULT);
  if (channel_id_group_id < 0) {
    fprintf(stderr, "Failed to open channel_id group %s\n", channel_id_group);
    exit(-1); // TODO(Haowen): fix this later
  }
  digitisation = GetFloatAttributeInGroup(channel_id_group_id, "digitisation");
  range = GetFloatAttributeInGroup(channel_id_group_id, "range");
  offset = GetFloatAttributeInGroup(channel_id_group_id, "offset");
  //scaling.sample_rate = GetFloatAttributeInGroup(channel_id_group_id, "sampling_rate");
  fast5_err = H5Gclose(channel_id_group_id);
  assert(fast5_err >= 0);
  // Get raw signal and convert it into current
  std::string read_signal_dataset = std::string(LEGACY_FAST5_RAW_ROOT) + "/" + std::string(group_name) + "/Signal";
  hid_t dataset_id = H5Dopen(fast5_file.hdf5_file, read_signal_dataset.c_str(), H5P_DEFAULT);
  if (dataset_id < 0) {
    fprintf(stderr, "Failed to open dataset '%s' to read signal from.\n", read_signal_dataset.c_str());
    return;
  }
  // Get an identifier for a copy of the dataspace for a dataset.
  hid_t dataspace_id = H5Dget_space(dataset_id);
  if (dataspace_id < 0) {
    fprintf(stderr, "Failed to create copy of dataspace for signal %s.\n", read_signal_dataset.c_str());
    goto cleanup3;
  }
  hsize_t first_dimension_size;
  num_dims = H5Sget_simple_extent_dims(dataspace_id, &first_dimension_size, NULL);
  assert(num_dims == 1); // There should be 1 dimension and the size of that dimension is the signal length
  signal_length = first_dimension_size;
  signal_values = (float*)calloc(signal_length, sizeof(float));
  fast5_err = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, signal_values);
  if (fast5_err < 0) {
    fprintf(stderr, "Failed to load read signal data from dataset %s.\n", read_signal_dataset.c_str());
    free(signal_values);
    goto cleanup4;
  }
  // convert to pA
  scale = range / digitisation;
  for (size_t i = 0; i < signal_length; i++) {
    signal_values[i] = (signal_values[i] + offset) * scale;
  }
  signals_.emplace_back(Signal{read_name, digitisation, range, offset, signal_length, signal_values, NULL});
cleanup4:
  H5Sclose(dataspace_id);
cleanup3:
  H5Dclose(dataset_id);
  //std::cerr << "Read name: " << name << ", # signal points: " << signal_length << ".\n";
}

void SignalBatch::NormalizeSignalAt(size_t signal_index) {
  // Should use a linear algorithm like median of medians
  // But for now let us use sort
  std::vector<float> tmp_signal((signals_[signal_index]).signal_values, (signals_[signal_index]).signal_values + (signals_[signal_index]).signal_length);
  std::nth_element(tmp_signal.begin(), tmp_signal.begin() + tmp_signal.size() / 2, tmp_signal.end());
  float signal_median = tmp_signal[tmp_signal.size() / 2]; // This is a fake median, but should be okay for a quick implementation
  for (size_t i = 0; i < tmp_signal.size(); ++i) {
    tmp_signal[i] = std::abs(tmp_signal[i] - signal_median);
  }
  std::nth_element(tmp_signal.begin(), tmp_signal.begin() + tmp_signal.size() / 2, tmp_signal.end());
  float MAD = tmp_signal[tmp_signal.size() / 2]; // Again, fake MAD, ok for a quick implementation
  // Now we can normalize signal
  for (size_t i = 0; i < signals_[signal_index].signal_length; ++i) {
    ((signals_[signal_index]).signal_values)[i] = (((signals_[signal_index]).signal_values)[i] - signal_median) / MAD;
  }
}

void SignalBatch::ConvertSequencesToSignals(const SequenceBatch &sequence_batch, const PoreModel &pore_model, size_t num_sequences) {
  double real_start_time = GetRealTime();
  for (size_t sequence_index = 0; sequence_index < num_sequences; ++sequence_index) {
    size_t sequence_length = sequence_batch.GetSequenceLengthAt(sequence_index);
    float *signal_values = pore_model.GetLevelMeansAt(sequence_batch.GetSequenceAt(sequence_index), 0, sequence_length);
    float *negative_signal_values = pore_model.GetLevelMeansAt(sequence_batch.GetNegativeSequenceAt(sequence_index).data(), 0, sequence_length);
    char *name =  (char*)calloc(1 + sequence_batch.GetSequenceNameLengthAt(sequence_index), sizeof(char));
    strcpy(name, sequence_batch.GetSequenceNameAt(sequence_index));
    signals_.emplace_back(Signal{name, 0, 0, 0, sequence_length - pore_model.GetKmerSize() + 1, signal_values, negative_signal_values});
  }
  std::cerr << "Convert " << num_sequences << " sequences to signals in " << GetRealTime() - real_start_time << "s.\n";
}
} // namespace sigmap
