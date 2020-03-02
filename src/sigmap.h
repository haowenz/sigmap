#ifndef SIGMAP_H_
#define SIGMAP_H_

#include <string>

namespace sigmap {
class SigmapDriver {
 public:
  SigmapDriver() {}
  ~SigmapDriver() {}
  void ParseArgsAndRun(int argc, char *argv[]);
};

class Sigmap {
 public:
  Sigmap(){}
  ~Sigmap(){}
  Sigmap(const std::string &reference_file_path, const std::string &pore_model_file_path, const std::string &signal_directory) : reference_file_path_(reference_file_path), pore_model_file_path_(pore_model_file_path), signal_directory_(signal_directory) {}
  void Map();
 protected:
  std::string reference_file_path_;
  std::string pore_model_file_path_;
  std::string signal_directory_;
  std::string output_path_;
};
} // namespace sigmap

#endif // SIGMAP_H_
