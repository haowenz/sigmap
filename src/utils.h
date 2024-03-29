#ifndef UTILS_H_
#define UTILS_H_

#include <assert.h>
#include <dirent.h>
#include <float.h>
#include <hdf5.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "cwt.h"

//#define DEBUG
#define LEGACY_FAST5_RAW_ROOT "/Raw/Reads/"

namespace sigmap {
struct FAST5File {
  hid_t hdf5_file;
  bool is_multi_fast5;
};

inline static bool IsDirectory(const std::string &dir_path) {
  auto dir = opendir(dir_path.c_str());
  if (not dir) {
    return false;
  }
  closedir(dir);
  return true;
}

inline static std::vector<std::string> ListDirectory(
    const std::string &dir_path) {
  std::vector<std::string> res;
  DIR *dir;
  struct dirent *ent;
  dir = opendir(dir_path.c_str());
  if (not dir) {
    return res;
  }
  while ((ent = readdir(dir)) != nullptr) {
    res.push_back(ent->d_name);
  }
  closedir(dir);
  return res;
}

inline static double GetRealTime() {
  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp, &tzp);
  return tp.tv_sec + tp.tv_usec * 1e-6;
}

inline static double GetCPUTime() {
  struct rusage r;
  getrusage(RUSAGE_SELF, &r);
  return r.ru_utime.tv_sec + r.ru_stime.tv_sec +
         1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

inline static void ExitWithMessage(const std::string &message) {
  std::cerr << message << std::endl;
  exit(-1);
}

// For sequence manipulation
static constexpr uint8_t char_to_uint8_table_[256] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
static constexpr char uint8_to_char_table_[8] = {'A', 'C', 'G', 'T',
                                                 'N', 'N', 'N', 'N'};

inline static uint8_t CharToUint8(const char c) {
  return char_to_uint8_table_[(uint8_t)c];
}

inline static char Uint8ToChar(const uint8_t i) {
  return uint8_to_char_table_[i];
}

inline static uint64_t GenerateSeedFromSequence(const char *sequence,
                                                size_t sequence_length,
                                                uint32_t start_position,
                                                uint32_t seed_length) {
  uint64_t mask = (((uint64_t)1) << (2 * seed_length)) - 1;
  uint64_t seed = 0;
  for (uint32_t i = 0; i < seed_length; ++i) {
    if (start_position + i < sequence_length) {
      uint8_t current_base = CharToUint8(sequence[i + start_position]);
      if (current_base < 4) {                        // not an ambiguous base
        seed = ((seed << 2) | current_base) & mask;  // forward k-mer
      } else {
        seed = (seed << 2) & mask;  // N->A
      }
    } else {
      seed = (seed << 2) & mask;  // Pad A
    }
  }
  return seed;
}

// For FAST5 manipulation
inline static FAST5File OpenFAST5(const std::string &fast5_file_path) {
  FAST5File fast5_file;
  fast5_file.hdf5_file =
      H5Fopen(fast5_file_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (fast5_file.hdf5_file < 0) {
    fprintf(stderr, "could not open fast5 file: %s\n", fast5_file_path.c_str());
  }
  // Check for attribute that indicates whether it is single or multi-fast5
  // see: https://community.nanoporetech.com/posts/multi-fast5-format
  const std::string indicator_p1 = "/UniqueGlobalKey/";
  const std::string indicator_p2 = indicator_p1 + "tracking_id/";
  bool has_indicator =
      H5Lexists(fast5_file.hdf5_file, indicator_p1.c_str(), H5P_DEFAULT) &&
      H5Lexists(fast5_file.hdf5_file, indicator_p2.c_str(), H5P_DEFAULT);
  fast5_file.is_multi_fast5 = !has_indicator;
  return fast5_file;
}

inline static void CloseFAST5(FAST5File &fast5_file) {
  H5Fclose(fast5_file.hdf5_file);
}

inline static float GetFloatAttributeInGroup(hid_t group_id,
                                             const char *attribute_name) {
  // The group_id should be checked somewhere else!
  float attribute_value = NAN;
  hid_t attribute_id = H5Aopen(group_id, attribute_name, H5P_DEFAULT);
  if (attribute_id < 0) {
    fprintf(stderr, "Failed to open attribute '%s' for reading.",
            attribute_name);
    return attribute_value;
  }
  herr_t err = H5Aread(attribute_id, H5T_NATIVE_FLOAT, &attribute_value);
  if (err < 0) {
    fprintf(stderr, "error reading attribute %s\n", attribute_name);
    exit(EXIT_FAILURE);
  }
  H5Aclose(attribute_id);
  return attribute_value;
}

inline static void GetStringAttributeInGroup(hid_t group_id,
                                             const char *attribute_name,
                                             char **string_attribute) {
  // The group_id should be checked somewhere else!
  hid_t attribute_id = H5Aopen(group_id, attribute_name, H5P_DEFAULT);
  hid_t attribute_type_id;
  hid_t native_type_id;
  htri_t is_variable_string;
  if (attribute_id < 0) {
    fprintf(stderr, "Failed to open attribute '%s' for reading.",
            attribute_name);
    *string_attribute = NULL;
    goto close_attr;
  }
  // Get data type and check it is a fixed-length string
  attribute_type_id = H5Aget_type(attribute_id);
  if (attribute_type_id < 0) {
    fprintf(stderr, "failed to get attribute type %s\n", attribute_name);
    *string_attribute = NULL;
    goto close_type;
  }
  if (H5Tget_class(attribute_type_id) != H5T_STRING) {
    fprintf(stderr, "attribute %s is not a string\n", attribute_name);
    *string_attribute = NULL;
    goto close_type;
  }
  native_type_id = H5Tget_native_type(attribute_type_id, H5T_DIR_ASCEND);
  if (native_type_id < 0) {
    fprintf(stderr, "failed to get native type for %s\n", attribute_name);
    *string_attribute = NULL;
    goto close_native_type;
  }
  is_variable_string = H5Tis_variable_str(attribute_type_id);
  if (is_variable_string > 0) {
    // variable length string
    // char* buffer;
    // herr_t err = H5Aread(attribute_id, native_type_id, &string_attribute);
    herr_t err = H5Aread(attribute_id, native_type_id, string_attribute);
    if (err < 0) {
      fprintf(stderr, "error reading attribute %s\n", attribute_name);
      exit(EXIT_FAILURE);
    }
    // out = buffer;
    // free(buffer);
    // buffer = NULL;
  } else if (is_variable_string == 0) {
    // fixed length string
    // Get the storage size and allocate
    size_t storage_size = H5Aget_storage_size(attribute_id);
    ;
    // char* buffer;
    // buffer = (char*)calloc(storage_size + 1, sizeof(char));
    *string_attribute = (char *)calloc(storage_size + 1, sizeof(char));
    // finally read the attribute
    herr_t err = H5Aread(attribute_id, attribute_type_id, *string_attribute);
    if (err < 0) {
      fprintf(stderr, "error reading attribute %s\n", attribute_name);
      exit(EXIT_FAILURE);
    }
    // clean up
    // free(buffer);
  } else {
    fprintf(stderr,
            "error when determing whether it is a variable string for "
            "attribute %s\n",
            attribute_name);
    exit(EXIT_FAILURE);
  }
  // H5Aread(attribute_id, H5T_NATIVE_FLOAT, &attribute_value);
  // H5Aclose(attribute_id);
close_native_type:
  H5Tclose(native_type_id);
close_type:
  H5Tclose(attribute_type_id);
close_attr:
  H5Aclose(attribute_id);
  // close_group:
  //  H5Gclose(group);
}

// For signal processing
static inline void GetSignalMinAndMax(const float *signal_values,
                                      size_t signal_length, float &min,
                                      float &max) {
  min = signal_values[0];
  max = signal_values[0];
  for (size_t position = 0; position < signal_length; ++position) {
    if (signal_values[position] < min) {
      min = signal_values[position];
    }
    if (signal_values[position] > max) {
      max = signal_values[position];
    }
  }
}

static inline size_t GetSignalsTotalLength(
    const std::vector<std::vector<float> > &signals) {
  size_t length = 0;
  for (size_t i = 0; i < signals.size(); ++i) {
    length += signals[i].size();
  }
  return length;
}

static inline void SaveVectorToFile(const float *signal_values,
                                    size_t signal_length,
                                    const std::string &file_path) {
  FILE *output_file = fopen(file_path.c_str(), "w");
  assert(output_file != NULL);
  for (size_t signal_position = 0; signal_position < signal_length - 1;
       ++signal_position) {
    fprintf(output_file, "%f\t", signal_values[signal_position]);
  }
  fprintf(output_file, "%f\n", signal_values[signal_length - 1]);
  fclose(output_file);
}

static inline void SaveVectorToFile(const size_t *positions,
                                    size_t signal_length,
                                    const std::string &file_path) {
  FILE *output_file = fopen(file_path.c_str(), "w");
  assert(output_file != NULL);
  for (size_t signal_position = 0; signal_position < signal_length - 1;
       ++signal_position) {
    fprintf(output_file, "%lu\t", positions[signal_position]);
  }
  fprintf(output_file, "%lu\n", positions[signal_length - 1]);
  fclose(output_file);
}
}  // namespace sigmap

#endif  // UTILS_H_
