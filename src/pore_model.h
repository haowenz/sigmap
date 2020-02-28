#ifndef POREMODEL_H_
#define POREMODEL_H_

#include <string>
#include <vector>

#include "utils.h"

namespace sigmap {
struct PoreModelParameters {
  uint16_t kmer;	
  float level_mean;
  float level_stdv;
  float sd_mean;
  float sd_stdv;	
  //float weight; // No need to store weight
};

class PoreModel {
 public:
  PoreModel(){}
  ~PoreModel(){}
  void Load(const std::string &pore_model_file_path);
  void Print();

 protected:
  int kmer_size_;
  std::vector<PoreModelParameters> pore_models_;
};
} // namespace sigmap

#endif // POREMODEL_H_
