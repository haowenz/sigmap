#ifndef SIGNALBATCH_H_
#define SIGNALBATCH_H_

#include <string>
#include <vector>

#include "hdf5_tools.hpp"
#include "pore_model.h"
#include "sequence_batch.h"
#include "utils.h"

namespace sigmap {
struct Signal {
  // read group name
  std::string id;
  // channel_id
  float digitisation;
  float range;
  float offset;
  // signal
  std::vector<float> signal_values;
  std::vector<float> negative_signal_values;
  size_t GetSignalLength() const {
    return signal_values.size();
  }
};

class SignalBatch {
 public:
  SignalBatch(){}
  ~SignalBatch(){}
  void InitializeLoading(const std::string &signal_directory);
  void FinalizeLoading();
  size_t LoadAllReadSignals();
  void AddSignal(const hdf5_tools::File &file, const std::string &raw_path, const std::string &ch_path);
  void AddSignalsFromFAST5(const std::string &fast5_file_path);
  void NormalizeSignalAt(size_t signal_index);
  void ConvertSequencesToSignals(const SequenceBatch &sequence_batch, const PoreModel &pore_model, size_t num_sequences);
  const Signal& GetSignalAt(size_t signal_index) const {
    return signals_[signal_index];
  }
  const char* GetSignalNameAt(size_t signal_index) const {
    return signals_[signal_index].id.data();
  }
  size_t GetSignalLengthAt(size_t signal_index) const {
    return signals_[signal_index].signal_values.size();
  }

 protected:
  std::string signal_directory_;
  std::vector<Signal> signals_;
};
} // namespace sigmap

#endif // SIGNALBATCH_H_
