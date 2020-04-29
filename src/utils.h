#ifndef UTILS_H_
#define UTILS_H_

#include <algorithm>
#include <assert.h>
#include <dirent.h>
#include <hdf5.h>
#include <iostream>
#include <cmath>
#include <float.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <vector>

#include "cwt.h"

#define LEGACY_FAST5_RAW_ROOT "/Raw/Reads/"

namespace sigmap {
struct FAST5File {
  hid_t hdf5_file;
  bool is_multi_fast5;
};

inline static bool IsDirectory(const std::string& dir_path) {
  auto dir = opendir(dir_path.c_str());
  if(not dir) {
    return false;
  }
  closedir(dir);
  return true;
}

inline static std::vector<std::string> ListDirectory(const std::string& dir_path) {
  std::vector<std::string> res;
  DIR* dir;
  struct dirent *ent;
  dir = opendir(dir_path.c_str());
  if(not dir) {
    return res;
  }
  while((ent = readdir(dir)) != nullptr) {
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
  return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

inline static void ExitWithMessage(const std::string &message) {
  std::cerr << message << std::endl;
  exit(-1);
}

// For sequence manipulation
static constexpr uint8_t char_to_uint8_table_[256] = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
static constexpr char uint8_to_char_table_[8] = {'A', 'C', 'G', 'T', 'N', 'N', 'N', 'N'};

inline static uint8_t CharToUint8(const char c) {
  return char_to_uint8_table_[(uint8_t)c];
}

inline static char Uint8ToChar(const uint8_t i) {
  return uint8_to_char_table_[i];
}

inline static uint64_t GenerateSeedFromSequence(const char *sequence, size_t sequence_length, uint32_t start_position, uint32_t seed_length) {
  //const char *sequence = GetSequenceAt(sequence_index);
  //uint32_t sequence_length = GetSequenceLengthAt(sequence_index);
  uint64_t mask = (((uint64_t)1) << (2 * seed_length)) - 1;
  uint64_t seed = 0;
  for (uint32_t i = 0; i < seed_length; ++i) {
    if (start_position + i < sequence_length) {
      uint8_t current_base = CharToUint8(sequence[i + start_position]);
      if (current_base < 4) { // not an ambiguous base
        seed = ((seed << 2) | current_base) & mask; // forward k-mer
      } else {
        seed = (seed << 2) & mask; // N->A
      }
    } else {
      seed = (seed << 2) & mask; // Pad A
    }
  }
  return seed;
}

// For FAST5 manipulation
inline static FAST5File OpenFAST5(const std::string &fast5_file_path) {
  FAST5File fast5_file;
  fast5_file.hdf5_file = H5Fopen(fast5_file_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (fast5_file.hdf5_file < 0) {
    fprintf(stderr, "could not open fast5 file: %s\n", fast5_file_path.c_str());
  }
  // Check for attribute that indicates whether it is single or multi-fast5
  // see: https://community.nanoporetech.com/posts/multi-fast5-format
  const std::string indicator_p1 = "/UniqueGlobalKey/";
  const std::string indicator_p2 = indicator_p1 + "tracking_id/";
  bool has_indicator = H5Lexists(fast5_file.hdf5_file, indicator_p1.c_str(), H5P_DEFAULT) && H5Lexists(fast5_file.hdf5_file, indicator_p2.c_str(), H5P_DEFAULT);
  fast5_file.is_multi_fast5 = !has_indicator;
  return fast5_file;
}

inline static void CloseFAST5(FAST5File& fast5_file) {
  H5Fclose(fast5_file.hdf5_file);
}

inline static float GetFloatAttributeInGroup(hid_t group_id, const char *attribute_name) {
  // The group_id should be checked somewhere else!
  float attribute_value = NAN;
  hid_t attribute_id = H5Aopen(group_id, attribute_name, H5P_DEFAULT);
  if (attribute_id < 0) {
    fprintf(stderr, "Failed to open attribute '%s' for reading.", attribute_name);
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

inline static void GetStringAttributeInGroup(hid_t group_id, const char *attribute_name, char **string_attribute) {
  // The group_id should be checked somewhere else!
  hid_t attribute_id = H5Aopen(group_id, attribute_name, H5P_DEFAULT);
  hid_t attribute_type_id;
  hid_t native_type_id;
  htri_t is_variable_string; 
  if (attribute_id < 0) {
    fprintf(stderr, "Failed to open attribute '%s' for reading.", attribute_name);
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
    //char* buffer;
    //herr_t err = H5Aread(attribute_id, native_type_id, &string_attribute);
    herr_t err = H5Aread(attribute_id, native_type_id, string_attribute);
    if (err < 0) {
      fprintf(stderr, "error reading attribute %s\n", attribute_name);
      exit(EXIT_FAILURE);
    }
    //out = buffer;
    //free(buffer);
    //buffer = NULL;
  } else if (is_variable_string == 0) {
    // fixed length string
    // Get the storage size and allocate
    size_t storage_size = H5Aget_storage_size(attribute_id);;
    //char* buffer;
    //buffer = (char*)calloc(storage_size + 1, sizeof(char));
    *string_attribute = (char*)calloc(storage_size + 1, sizeof(char));
    // finally read the attribute
    herr_t err = H5Aread(attribute_id, attribute_type_id, *string_attribute);
    if (err < 0) {
      fprintf(stderr, "error reading attribute %s\n", attribute_name);
      exit(EXIT_FAILURE);
    }
    // clean up
    //free(buffer);
  } else {
    fprintf(stderr, "error when determing whether it is a variable string for attribute %s\n", attribute_name);
    exit(EXIT_FAILURE);
  }
  //H5Aread(attribute_id, H5T_NATIVE_FLOAT, &attribute_value);
  //H5Aclose(attribute_id);
close_native_type:
  H5Tclose(native_type_id);    
close_type:
  H5Tclose(attribute_type_id);
close_attr:
  H5Aclose(attribute_id);
//close_group:
//  H5Gclose(group);
}

// For signal processing
static float GetNormalizedSignal(const float *signal, size_t signal_length, std::vector<float> &normalized_signal) {
  // Should use a linear algorithm like median of medians
  // One such better algorithm can be found here: https://rcoh.me/posts/linear-time-median-finding/
  // But for now let us use sort
  normalized_signal.assign(signal, signal + signal_length);
  std::nth_element(normalized_signal.begin(), normalized_signal.begin() + signal_length / 2, normalized_signal.end());
  float signal_median = normalized_signal[signal_length / 2]; // This is a fake median, but should be okay for a quick implementation
  for (size_t i = 0; i < signal_length; ++i) {
    normalized_signal[i] = std::abs(normalized_signal[i] - signal_median);
  }
  std::nth_element(normalized_signal.begin(), normalized_signal.begin() + signal_length / 2, normalized_signal.end());
  float MAD = normalized_signal[signal_length / 2]; // Again, fake MAD, ok for a quick implementation
  // Now we can normalize signal
  for (size_t i = 0; i < signal_length; ++i) {
    normalized_signal[i] = (signal[i] - signal_median) / MAD;
  }
  return MAD;
}

static float GetZscoreNormalizedSignal(const float *signal, size_t signal_length, std::vector<float> &normalized_signal) {
  // Calculate mean
  float mean = 0;
  for (size_t i = 0; i < signal_length; ++i) {
    mean += signal[i];
  }
  mean /= signal_length;
  // Calculate standard deviation
  float SD = 0;
  for (size_t i = 0; i < signal_length; ++i) {
    SD += (signal[i] - mean) * (signal[i] - mean);
  }
  SD /= (signal_length - 1);
  SD = sqrt(SD); 
  // Now we can normalize signal
  for (size_t i = 0; i < signal_length; ++i) {
    normalized_signal.emplace_back((signal[i] - mean) / SD);
  }
  return SD;
}

static void GetSignalMinAndMax(const float *signal, size_t signal_length, float &min, float &max) {
  min = signal[0];
  max = signal[0];
  for (size_t position = 0; position < signal_length; ++position) {
    if (signal[position] < min) {
      min = signal[position];
    }
    if (signal[position] > max) {
      max = signal[position];
    }
  }
}

static void GetCWTSignal(const float *signal, size_t signal_length, float scale0, std::vector<float> &cwt_signal) {
  char wave[] = "dog";
  char type[] = "pow";
  double param = 2.0;
  double dt = 1;
  double dj = 1; // Separation bewteen scales.
  int J = 1;
  int N = signal_length;
  cwt_object wt = cwt_init(wave, param, N, dt, J);
  setCWTScales(wt, scale0, dj, type, 2.0);
  cwt(wt, signal);
  for (size_t i = 0; i < signal_length; ++i) {
    cwt_signal.push_back(wt->output[i].re); 
  }
  //cwt_summary(wt);
  cwt_free(wt);
}

static void GetPeaks(const float *signal, size_t signal_length, float selective, std::vector<float> &peaks, std::vector<size_t> &peak_positions) {
  float previous_valley = signal[0];
  float previous_peak = signal[0];
  for (size_t i = 1; i < signal_length - 1; ++i) {
    if (signal[i] > signal[i - 1] && signal[i] >= signal[i + 1] && signal[i] >= previous_valley + selective) {
      peaks.push_back(signal[i]);
      peak_positions.push_back(i);
      previous_peak = signal[i];
    } else if (signal[i] < signal[i - 1] && signal[i] <= signal[i + 1] && signal[i] <= previous_peak - selective) {
      peaks.push_back(signal[i]);
      peak_positions.push_back(i);
      previous_valley = signal[i];
    }
  }
}

static size_t GetSignalsTotalLength(const std::vector<std::vector<float> > &signals) {
  size_t length = 0;
  for (size_t i = 0; i < signals.size(); ++i) {
    length += signals[i].size();
  }
  return length;
}

static void SaveVectorToFile(const float *signal, size_t signal_length, const std::string &file_path) {
  FILE *output_file = fopen(file_path.c_str(), "w");
  assert(output_file != NULL);
  for (size_t signal_position = 0; signal_position < signal_length - 1; ++signal_position) {
    fprintf(output_file, "%f\t", signal[signal_position]);
  }
  fprintf(output_file, "%f\n", signal[signal_length - 1]);
  fclose(output_file);
}

static void SaveVectorToFile(const size_t *positions, size_t signal_length, const std::string &file_path) {
  FILE *output_file = fopen(file_path.c_str(), "w");
  assert(output_file != NULL);
  for (size_t signal_position = 0; signal_position < signal_length - 1; ++signal_position) {
    fprintf(output_file, "%lu\t", positions[signal_position]);
  }
  fprintf(output_file, "%lu\n", positions[signal_length - 1]);
  fclose(output_file);
}

// scrappie code to do segmentation
struct Event {
  uint64_t start;
  size_t length;
  float mean;
  float stdv;
  //int pos;
  //int state;
};

//struct EventTable {
//  size_t n;
//  size_t start;
//  size_t end;
//  Event *event;
//};

//struct RawTable {
//  size_t n;
//  size_t start;
//  size_t end;
//  float *raw;
//};

struct DetectorArgs {
  size_t window_length1;
  size_t window_length2;
  float threshold1;
  float threshold2;
  float peak_height;
};

static DetectorArgs const event_detection_defaults = {
  .window_length1 = 3,
  .window_length2 = 6,
  .threshold1 = 1.4f,
  .threshold2 = 9.0f,
  .peak_height = 0.2f
};

static DetectorArgs const event_detection_rna = {
  .window_length1 = 7,
  .window_length2 = 14,
  .threshold1 = 2.5f,
  .threshold2 = 9.0f,
  .peak_height = 1.0f
};

struct Detector {
  int DEF_PEAK_POS;
  float DEF_PEAK_VAL;
  float *signal;
  size_t signal_length;
  float threshold;
  size_t window_length;
  size_t masked_to;
  int peak_pos;
  float peak_value;
  bool valid_peak;
};

/**
 *   Compute cumulative sum and sum of squares for a vector of data
 *
 *   Element i  sum (sumsq) is the sum (sum of squares) up to but
 *   excluding element i of the inputy data.
 *
 *   @param data      double[d_length]   Data to be summed over (in)
 *   @param sum       double[d_length + 1]   Vector to store sum (out)
 *   @param sumsq     double[d_length + 1]   Vector to store sum of squares (out)
 *   @param d_length                     Length of data vector
 **/
static inline void compute_sum_sumsq(const float *data, size_t d_length, std::vector<float> &sum, std::vector<float> &sumsq) {
  //RETURN_NULL_IF(NULL == data, );
  //RETURN_NULL_IF(NULL == sum, );
  //RETURN_NULL_IF(NULL == sumsq, );
  assert(d_length > 0);
  //sum[0] = 0.0f;
  //sumsq[0] = 0.0f;
  sum.emplace_back(0.0f);
  sumsq.emplace_back(0.0f);
  for (size_t i = 0; i < d_length; ++i) {
    //sum[i + 1] = sum[i] + data[i];
    //sumsq[i + 1] = sumsq[i] + data[i] * data[i];
    sum.emplace_back(sum[i] + data[i]);
    sumsq.emplace_back(sumsq[i] + data[i] * data[i]);
  }
}

/**
 *   Compute windowed t-statistic from summary information
 *
 *   @param sum       double[d_length]  Cumulative sums of data (in)
 *   @param sumsq     double[d_length]  Cumulative sum of squares of data (in)
 *   @param d_length                    Length of data vector
 *   @param w_length                    Window length to calculate t-statistic over
 *
 *   @returns float array containing tstats.  Returns NULL on error
 **/
static inline void compute_tstat(const float *sum, const float *sumsq, size_t d_length, size_t w_length, std::vector<float> &tstat) {
  assert(d_length > 0);
  assert(w_length > 0);
  //RETURN_NULL_IF(NULL == sum, NULL);
  //RETURN_NULL_IF(NULL == sumsq, NULL);

  //float *tstat = calloc(d_length, sizeof(float));
  tstat.resize(d_length);
  //RETURN_NULL_IF(NULL == tstat, NULL);

  const float eta = FLT_MIN;
  const float w_lengthf = (float)w_length;

  // Quick return:
  //   t-test not defined for number of points less than 2
  //   need at least as many points as twice the window length
  if (d_length < 2 * w_length || w_length < 2) {
    for (size_t i = 0; i < d_length; ++i) {
      //tstat[i] = 0.0f;
      tstat.emplace_back(0.0f);
    }
    return;
  }
  // fudge boundaries
  for (size_t i = 0; i < w_length; ++i) {
    tstat[i] = 0;
    tstat[d_length - i - 1] = 0;
  }

  // get to work on the rest
  for (size_t i = w_length; i <= d_length - w_length; ++i) {
    float sum1 = sum[i];
    float sumsq1 = sumsq[i];
    if (i > w_length) {
      sum1 -= sum[i - w_length];
      sumsq1 -= sumsq[i - w_length];
    }
    float sum2 = (float)(sum[i + w_length] - sum[i]);
    float sumsq2 = (float)(sumsq[i + w_length] - sumsq[i]);
    float mean1 = sum1 / w_lengthf;
    float mean2 = sum2 / w_lengthf;
    float combined_var = sumsq1 / w_lengthf - mean1 * mean1 + sumsq2 / w_lengthf - mean2 * mean2;

    // Prevent problem due to very small variances
    combined_var = fmaxf(combined_var, eta);

    //t-stat
    //  Formula is a simplified version of Student's t-statistic for the
    //  special case where there are two samples of equal size with
    //  differing variance
    const float delta_mean = mean2 - mean1;
    tstat[i] = fabs(delta_mean) / sqrt(combined_var / w_lengthf);
  }
}

/**
 *
 *   @returns array of length nsample whose elements contain peak positions
 *   Remaining elements are padded by zeros.
 **/
static inline void short_long_peak_detector(Detector *short_detector, Detector *long_detector, const float peak_height, std::vector<size_t> &peaks) {
  assert(short_detector->signal_length == long_detector->signal_length);
  //RETURN_NULL_IF(NULL == short_detector->signal, NULL);
  //RETURN_NULL_IF(NULL == long_detector->signal, NULL);

  const size_t ndetector = 2;
  Detector *detectors[ndetector] = {short_detector, long_detector};

  peaks.assign(short_detector->signal_length, 0);
  //size_t *peaks = calloc(short_detector->signal_length, sizeof(size_t));
  //RETURN_NULL_IF(NULL == peaks, NULL);

  size_t peak_count = 0;
  for (size_t i = 0; i < short_detector->signal_length; i++) {
    for (size_t k = 0; k < ndetector; k++) {
      Detector *detector = detectors[k];
      //Carry on if we've been masked out
      if (detector->masked_to >= i) {
        continue;
      }

      float current_value = detector->signal[i];

      if (detector->peak_pos == detector->DEF_PEAK_POS) {
        //CASE 1: We've not yet recorded a maximum
        if (current_value < detector->peak_value) {
          //Either record a deeper minimum...
          detector->peak_value = current_value;
        } else if (current_value - detector->peak_value > peak_height) {
          // ...or we've seen a qualifying maximum
          detector->peak_value = current_value;
          detector->peak_pos = i;
          //otherwise, wait to rise high enough to be considered a peak
        }
      } else {
        //CASE 2: In an existing peak, waiting to see if it is good
        if (current_value > detector->peak_value) {
          //Update the peak
          detector->peak_value = current_value;
          detector->peak_pos = i;
        }
        //Dominate other tstat signals if we're going to fire at some point
        if (detector == short_detector) {
          if (detector->peak_value > detector->threshold) {
            long_detector->masked_to = detector->peak_pos + detector->window_length;
            long_detector->peak_pos = long_detector->DEF_PEAK_POS;
            long_detector->peak_value = long_detector->DEF_PEAK_VAL;
            long_detector->valid_peak = false;
          }
        }
        //Have we convinced ourselves we've seen a peak
        if (detector->peak_value - current_value > peak_height && detector->peak_value > detector->threshold) {
          detector->valid_peak = true;
        }
        //Finally, check the distance if this is a good peak
        if (detector->valid_peak && (i - detector->peak_pos) > detector->window_length / 2) {
          //Emit the boundary and reset
          //peaks.emplace_back(detector->peak_pos);
          peaks[peak_count] = detector->peak_pos;
          peak_count++;
          detector->peak_pos = detector->DEF_PEAK_POS;
          detector->peak_value = current_value;
          detector->valid_peak = false;
        }
      }
    }
  }
  //return peaks;
}

/**  Create an event given boundaries
 *
 *   Note: Bounds are CADLAG (i.e. lower bound is contained in the interval but
 *   the upper bound is not).
 *
 *  @param start Index of lower bound
 *  @param end Index of upper bound
 *  @param sums
 *  @param sumsqs
 *  @param nsample  Total number of samples in read
 *
 *  @returns An initialised event.  A 'null' event is returned on error.
 **/
static inline Event create_event(size_t start, size_t end, const float *sums, const float *sumsqs, size_t nsample) {
  assert(start < nsample);
  assert(end <= nsample);

  Event event = { 0 };
  //event.pos = -1;
  //event.state = -1;
  //RETURN_NULL_IF(NULL == sums, event);
  //RETURN_NULL_IF(NULL == sumsqs, event);

  event.start = (uint64_t)start;
  event.length = (float)(end - start);
  event.mean = (float)(sums[end] - sums[start]) / event.length;
  float deltasqr = (sumsqs[end] - sumsqs[start]);
  float var = deltasqr / event.length - event.mean * event.mean;
  event.stdv = sqrtf(fmaxf(var, 0.0f));

  return event;
}

static inline void create_events(const size_t *peaks, const float *sums, const float *sumsqs, size_t nsample, std::vector<Event> &events) {
  //EventTable et = { 0 };
  //RETURN_NULL_IF(NULL == sums, et);
  //RETURN_NULL_IF(NULL == sumsqs, et);
  //RETURN_NULL_IF(NULL == peaks, et);

  // Count number of events found
  size_t n = 1;
  for (size_t i = 1; i < nsample; ++i) {
    if (peaks[i] > 0 && peaks[i] < nsample) {
      n++;
    }
  }

  //et.event = calloc(n, sizeof(Event));
  ////RETURN_NULL_IF(NULL == et.event, et);

  //et.n = n;
  //et.end = et.n;


  // First event -- starts at zero
  //et.event[0] = create_event(0, peaks[0], sums, sumsqs, nsample);
  events.emplace_back(create_event(0, peaks[0], sums, sumsqs, nsample));
  // Other events -- peak[i-1] -> peak[i]
  for (size_t ev = 1 ; ev < n - 1 ; ev++) {
    //et.event[ev] = create_event(peaks[ev - 1], peaks[ev], sums, sumsqs, nsample);
    events.emplace_back(create_event(peaks[ev - 1], peaks[ev], sums, sumsqs, nsample));
  }
  // Last event -- ends at nsample
  //et.event[n - 1] = create_event(peaks[n - 2], nsample, sums, sumsqs, nsample);
  events.emplace_back(create_event(peaks[n - 2], nsample, sums, sumsqs, nsample));

  //return et;
}

static inline void detect_events(const float *signal, size_t signal_length, const DetectorArgs &edparam, std::vector<Event> &events) {

  //EventTable et = { 0 };
  //RETURN_NULL_IF(NULL == rt.raw, et);

  std::vector<float> sums;// = calloc(signal.signal_length + 1, sizeof(double));
  sums.reserve(signal_length + 1);
  std::vector<float> sumsqs;// = calloc(signal_length + 1, sizeof(double));
  sumsqs.reserve(signal_length + 1);

  compute_sum_sumsq(signal, signal_length, sums, sumsqs);
  std::vector<float> tstat1;
  std::vector<float> tstat2;
  compute_tstat(sums.data(), sumsqs.data(), signal_length, edparam.window_length1, tstat1);
  compute_tstat(sums.data(), sumsqs.data(), signal_length, edparam.window_length2, tstat2);
  //float *tstat1 = compute_tstat(sums, sumsqs, rt.n, edparam.window_length1);
  //float *tstat2 = compute_tstat(sums, sumsqs, rt.n, edparam.window_length2);

  Detector short_detector = {
    .DEF_PEAK_POS = -1,
    .DEF_PEAK_VAL = FLT_MAX,
    .signal = tstat1.data(),
    .signal_length = signal_length,
    .threshold = edparam.threshold1,
    .window_length = edparam.window_length1,
    .masked_to = 0,
    .peak_pos = -1,
    .peak_value = FLT_MAX,
    .valid_peak = false
  };

  Detector long_detector = {
    .DEF_PEAK_POS = -1,
    .DEF_PEAK_VAL = FLT_MAX,
    .signal = tstat2.data(),
    .signal_length = signal_length,
    .threshold = edparam.threshold2,
    .window_length = edparam.window_length2,
    .masked_to = 0,
    .peak_pos = -1,
    .peak_value = FLT_MAX,
    .valid_peak = false
  };

  //size_t *peaks = short_long_peak_detector(&short_detector, &long_detector, edparam.peak_height);
  std::vector<size_t> peaks;
  short_long_peak_detector(&short_detector, &long_detector, edparam.peak_height, peaks);

  create_events(peaks.data(), sums.data(), sumsqs.data(), signal_length, events);
  std::cerr << "Detected " << events.size() << " events.\n";
  //et = create_events(peaks, sums, sumsqs, rt.n);

  //free(peaks);
  //free(tstat2);
  //free(tstat1);
  //free(sumsqs);
  //free(sums);

  //return et;
}
} // namespace sigmap

#endif // UTILS_H_
