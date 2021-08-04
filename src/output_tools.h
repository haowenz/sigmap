#ifndef OUTPUTTOOLS_H_
#define OUTPUTTOOLS_H_

#include <assert.h>

#include <cstring>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "sequence_batch.h"

namespace sigmap {

struct PAFMapping {
  uint32_t read_id;
  std::string read_name;
  uint32_t read_length;
  uint32_t read_start_position;
  uint32_t read_end_position;
  uint32_t fragment_start_position;
  uint32_t fragment_length;
  uint8_t mapq : 6, direction : 1, is_unique : 1;
  std::string tags;
  bool operator<(const PAFMapping &m) const {
    return std::tie(fragment_start_position, fragment_length, mapq, direction,
                    is_unique, read_id, read_start_position,
                    read_end_position) <
           std::tie(m.fragment_start_position, m.fragment_length, m.mapq,
                    m.direction, m.is_unique, m.read_id, m.read_start_position,
                    m.read_end_position);
  }
  bool operator==(const PAFMapping &m) const {
    return std::tie(fragment_start_position) ==
           std::tie(m.fragment_start_position);
  }
};

struct MappingWithBarcode {
  uint32_t read_id;
  uint32_t cell_barcode;
  uint32_t fragment_start_position;
  uint16_t fragment_length;
  uint8_t mapq : 6, direction : 1, is_unique : 1;
  bool operator<(const MappingWithBarcode &m) const {
    return std::tie(fragment_start_position, fragment_length, cell_barcode,
                    mapq, direction, is_unique, read_id) <
           std::tie(m.fragment_start_position, m.fragment_length,
                    m.cell_barcode, m.mapq, m.direction, m.is_unique,
                    m.read_id);
  }
  bool operator==(const MappingWithBarcode &m) const {
    return std::tie(cell_barcode, fragment_start_position) ==
           std::tie(m.cell_barcode, m.fragment_start_position);
  }
};

struct MappingWithoutBarcode {
  uint32_t read_id;
  uint32_t fragment_start_position;
  uint16_t fragment_length;
  uint8_t mapq : 6, direction : 1, is_unique : 1;
  bool operator<(const MappingWithoutBarcode &m) const {
    return std::tie(fragment_start_position, fragment_length, mapq, direction,
                    is_unique, read_id) <
           std::tie(m.fragment_start_position, m.fragment_length, m.mapq,
                    m.direction, m.is_unique, m.read_id);
  }
  bool operator==(const MappingWithoutBarcode &m) const {
    return std::tie(fragment_start_position) ==
           std::tie(m.fragment_start_position);
  }
};

template <typename MappingRecord>
struct TempMappingFileHandle {
  std::string file_path;
  FILE *file;
  bool all_loaded;
  uint32_t num_mappings;
  uint32_t block_size;
  uint32_t current_rid;
  uint32_t current_mapping_index;
  uint32_t num_mappings_on_current_rid;
  uint32_t num_loaded_mappings_on_current_rid;
  std::vector<MappingRecord>
      mappings;  // this vector only keep mappings on the same ref seq
  inline void InitializeTempMappingLoading(uint32_t num_reference_sequences) {
    file = fopen(file_path.c_str(), "rb");
    assert(file != NULL);
    all_loaded = false;
    current_rid = 0;
    fread(&num_mappings_on_current_rid, sizeof(size_t), 1, file);
    mappings.resize(block_size);
    num_loaded_mappings_on_current_rid = 0;
    std::cerr << "Block size: " << block_size << ", initialize temp file "
              << file_path << "\n";
  }
  inline void FinalizeTempMappingLoading() { fclose(file); }
  inline void LoadTempMappingBlock(uint32_t num_reference_sequences) {
    num_mappings = 0;
    while (num_mappings == 0) {
      // Only keep mappings on one ref seq, which means # mappings in buffer can
      // be less than block size Two cases: current ref seq has remainings or
      // not
      if (num_loaded_mappings_on_current_rid < num_mappings_on_current_rid) {
        // Check if # remains larger than block size
        uint32_t num_mappings_to_load_on_current_rid =
            num_mappings_on_current_rid - num_loaded_mappings_on_current_rid;
        if (num_mappings_to_load_on_current_rid > block_size) {
          num_mappings_to_load_on_current_rid = block_size;
        }
        // std::cerr << num_mappings_to_load_on_current_rid << " " <<
        // num_loaded_mappings_on_current_rid << " " <<
        // num_mappings_on_current_rid << "\n"; std::cerr << mappings.size() <<
        // "\n";
        fread(mappings.data(), sizeof(MappingRecord),
              num_mappings_to_load_on_current_rid, file);
        // std::cerr << "Load mappings\n";
        num_loaded_mappings_on_current_rid +=
            num_mappings_to_load_on_current_rid;
        num_mappings = num_mappings_to_load_on_current_rid;
      } else {
        // Move to next rid
        ++current_rid;
        if (current_rid < num_reference_sequences) {
          // std::cerr << "Load size\n";
          fread(&num_mappings_on_current_rid, sizeof(size_t), 1, file);
          // std::cerr << "Load size " << num_mappings_on_current_rid << "\n";
          num_loaded_mappings_on_current_rid = 0;
        } else {
          all_loaded = true;
          break;
        }
      }
    }
    current_mapping_index = 0;
  }
  inline void Next(uint32_t num_reference_sequences) {
    ++current_mapping_index;
    if (current_mapping_index >= num_mappings) {
      LoadTempMappingBlock(num_reference_sequences);
    }
  }
};

template <typename MappingRecord>
class OutputTools {
 public:
  OutputTools() {}
  virtual ~OutputTools() {}
  inline void OutputTempMapping(
      const std::string &temp_mapping_output_file_path,
      uint32_t num_reference_sequences,
      const std::vector<std::vector<MappingRecord> > &mappings) {
    FILE *temp_mapping_output_file =
        fopen(temp_mapping_output_file_path.c_str(), "wb");
    assert(temp_mapping_output_file != NULL);
    for (size_t ri = 0; ri < num_reference_sequences; ++ri) {
      // make sure mappings[ri] exists even if its size is 0
      size_t num_mappings = mappings[ri].size();
      fwrite(&num_mappings, sizeof(size_t), 1, temp_mapping_output_file);
      if (mappings[ri].size() > 0) {
        fwrite(mappings[ri].data(), sizeof(MappingRecord), mappings[ri].size(),
               temp_mapping_output_file);
      }
    }
    fclose(temp_mapping_output_file);
  }
  inline void LoadBinaryTempMapping(
      const std::string &temp_mapping_file_path,
      uint32_t num_reference_sequences,
      std::vector<std::vector<MappingRecord> > &mappings) {
    FILE *temp_mapping_file = fopen(temp_mapping_file_path.c_str(), "rb");
    assert(temp_mapping_file != NULL);
    for (size_t ri = 0; ri < num_reference_sequences; ++ri) {
      size_t num_mappings = 0;
      fread(&num_mappings, sizeof(size_t), 1, temp_mapping_file);
      if (num_mappings > 0) {
        mappings.emplace_back(std::vector<MappingRecord>(num_mappings));
        fread(&(mappings[ri].data()), sizeof(MappingRecord), num_mappings,
              temp_mapping_file);
      } else {
        mappings.emplace_back(std::vector<MappingRecord>());
      }
    }
    fclose(temp_mapping_file);
  }
  inline void InitializeMappingOutput(
      const std::string &mapping_output_file_path) {
    mapping_output_file_path_ = mapping_output_file_path;
    mapping_output_file_ = fopen(mapping_output_file_path_.c_str(), "w");
    assert(mapping_output_file_ != NULL);
  }
  inline void FinalizeMappingOutput() { fclose(mapping_output_file_); }
  inline void AppendMappingOutput(const std::string &line) {
    fprintf(mapping_output_file_, "%s", line.data());
  }
  inline void AppendUnmappedRead(uint32_t rid, const SequenceBatch &reference,
                                 const PAFMapping &mapping) {
    // const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    // uint32_t reference_sequence_length = reference.GetSequenceLengthAt(rid);
    // uint32_t mapping_end_position = mapping.fragment_start_position +
    // mapping.fragment_length;
    this->AppendMappingOutput(
        mapping.read_name + "\t" + std::to_string(mapping.read_length) +
        "\t*\t*\t*\t*\t*\t*\t*\t*\t*\t" + std::to_string(mapping.mapq) + "\t" +
        mapping.tags + "\n");
  }
  virtual void AppendMapping(uint32_t rid, const SequenceBatch &reference,
                             const MappingRecord &mapping) = 0;
  inline std::string GeneratePAFLine(
      const SequenceBatch &query_batch, uint32_t query_index,
      const int query_start, const int query_end, const char relative_strand,
      const SequenceBatch &target_batch, uint32_t target_index,
      const int target_start, const int target_end, const int num_matches,
      const int alignment_length, const int mapping_quality) {
    return std::string(query_batch.GetSequenceNameAt(query_index)) + "\t" +
           std::to_string(query_batch.GetSequenceLengthAt(query_index)) + "\t" +
           std::to_string(query_start) + "\t" + std::to_string(query_end) +
           "\t" + relative_strand + "\t" +
           std::string(target_batch.GetSequenceNameAt(target_index)) + "\t" +
           std::to_string(target_batch.GetSequenceLengthAt(target_index)) +
           "\t" + std::to_string(target_start) + "\t" +
           std::to_string(target_end) + "\t" + std::to_string(num_matches) +
           "\t" + std::to_string(alignment_length) + "\t" +
           std::to_string(mapping_quality) + "\n";
  }
  inline uint32_t GetNumMappings() const { return num_mappings_; }
  inline std::string Seed2Sequence(uint64_t seed, uint32_t seed_length) const {
    std::string sequence;
    sequence.reserve(seed_length);
    uint64_t mask = 3;
    for (uint32_t i = 0; i < seed_length; ++i) {
      sequence.push_back(SequenceBatch::Uint8ToChar(
          (seed >> ((seed_length - 1 - i) * 2)) & mask));
    }
    return sequence;
  }

  // inline void InitializeMatrixOutput(const std::string &matrix_output_prefix)
  // {
  //  matrix_output_prefix_ = matrix_output_prefix;
  //  matrix_output_file_ = fopen((matrix_output_prefix_ +
  //  "_matrix.mtx").c_str(), "w"); assert(matrix_output_file_ != NULL);
  //  peak_output_file_ = fopen((matrix_output_prefix_ + "_peaks.bed").c_str(),
  //  "w"); assert(peak_output_file_ != NULL); barcode_output_file_ =
  //  fopen((matrix_output_prefix_ + "_barcode.tsv").c_str(), "w");
  //  assert(barcode_output_file_ != NULL);
  //}
  // void OutputPeaks(uint32_t bin_size, uint32_t num_sequences, const
  // SequenceBatch &reference) {
  //  for (uint32_t rid = 0; rid < num_sequences; ++rid) {
  //    uint32_t sequence_length = reference.GetSequenceLengthAt(rid);
  //    const char *sequence_name = reference.GetSequenceNameAt(rid);
  //    for (uint32_t position = 0; position < sequence_length; position +=
  //    bin_size) {
  //      fprintf(peak_output_file_, "%s\t%u\t%u\n", sequence_name, position +
  //      1, position + bin_size);
  //    }
  //  }
  //}
  // void OutputPeaks(uint32_t peak_start_position, uint16_t peak_length,
  // uint32_t rid, const SequenceBatch &reference) {
  //  const char *sequence_name = reference.GetSequenceNameAt(rid);
  //  fprintf(peak_output_file_, "%s\t%u\t%u\n", sequence_name,
  //  peak_start_position + 1, peak_start_position + peak_length);
  //}
  // void AppendBarcodeOutput(uint32_t barcode_key) {
  //  fprintf(barcode_output_file_, "%s-1\n", Seed2Sequence(barcode_key,
  //  cell_barcode_length_).data());
  //}
  // void WriteMatrixOutputHead(uint64_t num_peaks, uint64_t num_barcodes,
  // uint64_t num_lines) {
  //  fprintf(matrix_output_file_, "%lu\t%lu\t%lu\n", num_peaks, num_barcodes,
  //  num_lines);
  //}
  // void AppendMatrixOutput(uint32_t peak_index, uint32_t barcode_index,
  // uint32_t num_mappings) {
  //  fprintf(matrix_output_file_, "%u\t%u\t%u\n", peak_index, barcode_index,
  //  num_mappings);
  //}
  inline void FinalizeMatrixOutput() {
    fclose(matrix_output_file_);
    fclose(peak_output_file_);
    fclose(barcode_output_file_);
  }

 protected:
  std::string mapping_output_file_path_;
  FILE *mapping_output_file_;
  uint32_t num_mappings_;
  uint32_t cell_barcode_length_ = 16;
  std::string matrix_output_prefix_;
  FILE *peak_output_file_;
  FILE *barcode_output_file_;
  FILE *matrix_output_file_;
};

// template <typename MappingRecord>
// class BEDOutputTools : public OutputTools<MappingRecord> {
//  inline void AppendMapping(uint32_t rid, const SequenceBatch &reference,
//  const MappingRecord &mapping) {
//    std::string strand = (mapping.direction & 1) == 1 ? "+" : "-";
//    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
//    uint32_t mapping_end_position = mapping.fragment_start_position +
//    mapping.fragment_length;
//    this->AppendMappingOutput(std::string(reference_sequence_name) + "\t" +
//    std::to_string(mapping.fragment_start_position) + "\t" +
//    std::to_string(mapping_end_position) + "\tN\t1000\t" + strand + "\n");
//  }
//};
//
// template <>
// class BEDOutputTools<MappingWithBarcode> : public
// OutputTools<MappingWithBarcode> {
//  inline void AppendMapping(uint32_t rid, const SequenceBatch &reference,
//  const MappingWithBarcode &mapping) {
//    std::string strand = (mapping.direction & 1) == 1 ? "+" : "-";
//    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
//    uint32_t mapping_end_position = mapping.fragment_start_position +
//    mapping.fragment_length;
//    this->AppendMappingOutput(std::string(reference_sequence_name) + "\t" +
//    std::to_string(mapping.fragment_start_position) + "\t" +
//    std::to_string(mapping_end_position) + "\t" +
//    Seed2Sequence(mapping.cell_barcode, cell_barcode_length_) +"\n");
//  }
//};

template <typename MappingRecord>
class PAFOutputTools : public OutputTools<MappingRecord> {};

template <>
class PAFOutputTools<PAFMapping> : public OutputTools<PAFMapping> {
  inline void AppendMapping(uint32_t rid, const SequenceBatch &reference,
                            const PAFMapping &mapping) {
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t reference_sequence_length = reference.GetSequenceLengthAt(rid);
    std::string strand = (mapping.direction & 1) == 1 ? "+" : "-";
    uint32_t mapping_end_position =
        mapping.fragment_start_position + mapping.fragment_length;
    this->AppendMappingOutput(
        mapping.read_name + "\t" + std::to_string(mapping.read_length) + "\t" +
        std::to_string(mapping.read_start_position) + "\t" +
        std::to_string(mapping.read_end_position) + "\t" + strand + "\t" +
        std::string(reference_sequence_name) + "\t" +
        std::to_string(reference_sequence_length) + "\t" +
        std::to_string(mapping.fragment_start_position) + "\t" +
        std::to_string(mapping_end_position) + "\t" +
        std::to_string(mapping.read_length) + "\t" +
        std::to_string(mapping.fragment_length) + "\t" +
        std::to_string(mapping.mapq) + "\t" + mapping.tags + "\n");
  }
};

template <>
class PAFOutputTools<MappingWithBarcode>
    : public OutputTools<MappingWithBarcode> {
  inline void AppendMapping(uint32_t rid, const SequenceBatch &reference,
                            const MappingWithBarcode &mapping) {
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t reference_sequence_length = reference.GetSequenceLengthAt(rid);
    std::string strand = (mapping.direction & 1) == 1 ? "+" : "-";
    uint32_t mapping_end_position =
        mapping.fragment_start_position + mapping.fragment_length;
    this->AppendMappingOutput(
        std::to_string(mapping.read_id) + "\t" +
        std::to_string(mapping.fragment_length) + "\t" + std::to_string(0) +
        "\t" + std::to_string(mapping.fragment_length) + "\t" + strand + "\t" +
        std::string(reference_sequence_name) + "\t" +
        std::to_string(reference_sequence_length) + "\t" +
        std::to_string(mapping.fragment_start_position) + "\t" +
        std::to_string(mapping_end_position) + "\t" +
        std::to_string(mapping.fragment_length) + "\t" +
        std::to_string(mapping.fragment_length) + "\t" +
        std::to_string(mapping.mapq) + "\n");
  }
};

template <>
class PAFOutputTools<MappingWithoutBarcode>
    : public OutputTools<MappingWithoutBarcode> {
  inline void AppendMapping(uint32_t rid, const SequenceBatch &reference,
                            const MappingWithoutBarcode &mapping) {
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t reference_sequence_length = reference.GetSequenceLengthAt(rid);
    std::string strand = (mapping.direction & 1) == 1 ? "+" : "-";
    uint32_t mapping_end_position =
        mapping.fragment_start_position + mapping.fragment_length;
    this->AppendMappingOutput(
        std::to_string(mapping.read_id) + "\t" +
        std::to_string(mapping.fragment_length) + "\t" + std::to_string(0) +
        "\t" + std::to_string(mapping.fragment_length) + "\t" + strand + "\t" +
        std::string(reference_sequence_name) + "\t" +
        std::to_string(reference_sequence_length) + "\t" +
        std::to_string(mapping.fragment_start_position) + "\t" +
        std::to_string(mapping_end_position) + "\t" +
        std::to_string(mapping.fragment_length) + "\t" +
        std::to_string(mapping.fragment_length) + "\t" +
        std::to_string(mapping.mapq) + "\n");
  }
};
}  // namespace sigmap

#endif  // OUTPUTTOOLS_H_
