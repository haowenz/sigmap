#ifndef SIGNALBATCH_H_
#define SIGNALBATCH_H_

#include <string>
#include <vector>

#include "utils.h"

namespace sigmap {
struct Signal {
  // read group name
  size_t name_length;
  char *name;
  // channel_id
  int digitisation;
  float range;
  int offset;
  // signal
  size_t signal_length;
  float *signal;
  Signal() {
    name = nullptr;
    signal = nullptr;
  }
  ~Signal() {
    if (name != nullptr) {
      free(name);
    }
    if (signal != nullptr) {
      free(signal);
    }
  }
};

class SignalBatch {
 public:
  SignalBatch(){}
  ~SignalBatch(){}
  void InitializeLoading(const std::string &signal_directory);
  void FinalizeLoading();
  void LoadAllReadSignals();
  void AddSignalsFromFAST5(const std::string &fast5_file_path);
  void AddSignalFromSingleFAST5(const FAST5File& fast5_file);

 protected:
  std::string signal_directory_;
  std::vector<Signal> signals_;
};
} // namespace sigmap

#endif // SIGNALBATCH_H_
