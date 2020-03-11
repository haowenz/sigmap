#ifndef SIGNALBATCH_H_
#define SIGNALBATCH_H_

#include <string>
#include <vector>

#include "pore_model.h"
#include "sequence_batch.h"
#include "utils.h"

namespace sigmap {
struct Signal {
  // read group name
  char *name;
  // channel_id
  float digitisation;
  float range;
  float offset;
  // signal
  size_t signal_length;
  float *signal;
};

class SignalBatch {
 public:
  SignalBatch(){}
  ~SignalBatch() {
    for (size_t read_index = 0; read_index < signals_.size(); ++read_index) {
      if (signals_[read_index].name != nullptr) {
        free(signals_[read_index].name);
      }
      if (signals_[read_index].signal != nullptr) {
        free(signals_[read_index].signal);
      }
    }
  }
  void InitializeLoading(const std::string &signal_directory);
  void FinalizeLoading();
  size_t LoadAllReadSignals();
  void AddSignalsFromFAST5(const std::string &fast5_file_path);
  void AddSignalFromSingleFAST5(const FAST5File& fast5_file);
  void NormalizeSignalAt(size_t signal_index);
  void ConvertSequencesToSignals(const SequenceBatch &sequence_batch, const PoreModel &pore_model, size_t num_sequences);
  const Signal& GetSignalAt(size_t signal_index) const {
    return signals_[signal_index];
  }
  const char* GetSignalNameAt(size_t signal_index) const {
    return signals_[signal_index].name;
  }

 protected:
  std::string signal_directory_;
  std::vector<Signal> signals_;
};
} // namespace sigmap

#endif // SIGNALBATCH_H_
