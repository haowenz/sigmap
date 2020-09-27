#ifndef EVENT_H_
#define EVENT_H_

#include <assert.h>
#include <cmath>
#include <float.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <vector>

#include "utils.h"

namespace sigmap {
struct Event {
  uint64_t start;
  size_t length;
  float mean;
  float stdv;
};

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
  float *signal_values;
  size_t signal_length;
  float threshold;
  size_t window_length;
  size_t masked_to;
  int peak_pos;
  float peak_value;
  bool valid_peak;
};

static inline void ComputePrefixSumAndPrefixSumSquares(const float *data, size_t data_length, std::vector<float> &prefix_sum, std::vector<float> &prefix_sum_square) {
  assert(data_length > 0);
  prefix_sum.emplace_back(0.0f);
  prefix_sum_square.emplace_back(0.0f);
  for (size_t i = 0; i < data_length; ++i) {
    prefix_sum.emplace_back(prefix_sum[i] + data[i]);
    prefix_sum_square.emplace_back(prefix_sum_square[i] + data[i] * data[i]);
  }
}

static inline void ComputeTStat(const float *prefix_sum, const float *prefix_sum_square, size_t signal_length, size_t window_length, std::vector<float> &tstat) {
  const float eta = FLT_MIN;
  // Quick return:
  // t-test not defined for number of points less than 2
  // need at least as many points as twice the window length
  if (signal_length < 2 * window_length || window_length < 2) {
    for (size_t i = 0; i < signal_length; ++i) {
      tstat.emplace_back(0.0f);
    }
    return;
  }
  // fudge boundaries
  for (size_t i = 0; i < window_length; ++i) {
    tstat.emplace_back(0.0f);
  }
  // get to work on the rest
  for (size_t i = window_length; i <= signal_length - window_length; ++i) {
    float sum1 = prefix_sum[i];
    float sumsq1 = prefix_sum_square[i];
    if (i > window_length) {
      sum1 -= prefix_sum[i - window_length];
      sumsq1 -= prefix_sum_square[i - window_length];
    }
    float sum2 = prefix_sum[i + window_length] - prefix_sum[i];
    float sumsq2 = prefix_sum_square[i + window_length] - prefix_sum_square[i];
    float mean1 = sum1 / window_length;
    float mean2 = sum2 / window_length;
    float combined_var = sumsq1 / window_length - mean1 * mean1 + sumsq2 / window_length - mean2 * mean2;
    // Prevent problem due to very small variances
    combined_var = fmaxf(combined_var, eta);
    //t-stat
    //  Formula is a simplified version of Student's t-statistic for the
    //  special case where there are two samples of equal size with
    //  differing variance
    const float delta_mean = mean2 - mean1;
    tstat.emplace_back(fabs(delta_mean) / sqrt(combined_var / window_length));
  }
  // fudge boundaries
  for (size_t i = 0; i < window_length; ++i) {
    tstat.emplace_back(0.0f);
  }
}

static inline void GeneratePeaksUsingMultiWindows(Detector *short_detector, Detector *long_detector, const float peak_height, std::vector<size_t> &peaks) {
  assert(short_detector->signal_length == long_detector->signal_length);
  const size_t ndetector = 2;
  Detector *detectors[ndetector] = {short_detector, long_detector};
  peaks.reserve(short_detector->signal_length);
  for (size_t i = 0; i < short_detector->signal_length; i++) {
    for (size_t k = 0; k < ndetector; k++) {
      Detector *detector = detectors[k];
      //Carry on if we've been masked out
      if (detector->masked_to >= i) {
        continue;
      }
      float current_value = detector->signal_values[i];
      if (detector->peak_pos == detector->DEF_PEAK_POS) {
        //CASE 1: We've not yet recorded a maximum
        if (current_value < detector->peak_value) {
          //Either record a deeper minimum...
          detector->peak_value = current_value;
        } else if (current_value - detector->peak_value > peak_height) { // TODO(Haowen): this might cause overflow, need to fix this
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
          peaks.emplace_back(detector->peak_pos);
          detector->peak_pos = detector->DEF_PEAK_POS;
          detector->peak_value = current_value;
          detector->valid_peak = false;
        }
      }
    }
  }
}

static inline Event CreateEvent(size_t start, size_t end, const float *prefix_sum, const float *prefix_sum_square, size_t signal_length) {
  assert(start < signal_length);
  assert(end <= signal_length);
  Event event;
  event.start = start;
  event.length = end - start;
  event.mean = (prefix_sum[end] - prefix_sum[start]) / event.length;
  float deltasqr = prefix_sum_square[end] - prefix_sum_square[start];
  float var = deltasqr / event.length - event.mean * event.mean;
  event.stdv = sqrtf(fmaxf(var, 0.0f));
  return event;
}

static inline void CreateEvents(const size_t *peaks, uint32_t peak_size, const float *prefix_sum, const float *prefix_sum_square, size_t signal_length, std::vector<Event> &events) {
  // Count number of events found
  size_t num_events = 1;
  for (size_t i = 1; i < peak_size; ++i) {
    if (peaks[i] > 0 && peaks[i] < signal_length) {
      num_events++;
    }
  }
  // First event -- starts at zero
  events.emplace_back(CreateEvent(0, peaks[0], prefix_sum, prefix_sum_square, signal_length));
  // Other events -- peak[i-1] -> peak[i]
  for (size_t pi = 1 ; pi < num_events - 1 ; pi++) {
    events.emplace_back(CreateEvent(peaks[pi - 1], peaks[pi], prefix_sum, prefix_sum_square, signal_length));
  }
  // Last event -- ends at signal_length
  events.emplace_back(CreateEvent(peaks[num_events - 2], signal_length, prefix_sum, prefix_sum_square, signal_length));
}

static inline void DetectEvents(const float *signal_values, size_t signal_length, const DetectorArgs &edparam, std::vector<float> &prefix_sum, std::vector<float> &prefix_sum_square, std::vector<float> &tstat1, std::vector<float> &tstat2, std::vector<size_t> &peaks, std::vector<Event> &events) {
  prefix_sum.reserve(signal_length + 1);
  prefix_sum_square.reserve(signal_length + 1);
  ComputePrefixSumAndPrefixSumSquares(signal_values, signal_length, prefix_sum, prefix_sum_square);
  ComputeTStat(prefix_sum.data(), prefix_sum_square.data(), signal_length, edparam.window_length1, tstat1);
  ComputeTStat(prefix_sum.data(), prefix_sum_square.data(), signal_length, edparam.window_length2, tstat2);
  Detector short_detector = {
    .DEF_PEAK_POS = -1,
    .DEF_PEAK_VAL = FLT_MAX,
    .signal_values = tstat1.data(),
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
    .signal_values = tstat2.data(),
    .signal_length = signal_length,
    .threshold = edparam.threshold2,
    .window_length = edparam.window_length2,
    .masked_to = 0,
    .peak_pos = -1,
    .peak_value = FLT_MAX,
    .valid_peak = false
  };
  GeneratePeaksUsingMultiWindows(&short_detector, &long_detector, edparam.peak_height, peaks);
  CreateEvents(peaks.data(), peaks.size(), prefix_sum.data(), prefix_sum_square.data(), signal_length, events);
#ifdef DEBUG
  std::cerr << "Detected " << events.size() << " events.\n";
#endif
}
} // namespace sigmap

#endif // EVENT_H_
