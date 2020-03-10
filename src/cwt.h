#ifndef CWT_H_
#define CWT_H_

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>

namespace sigmap {
// FFT
#ifndef fft_type
#define fft_type float
#endif

#define PI2 6.28318530717958647692528676655900577

typedef struct fft_set* fft_object;

typedef struct fft_t {
  fft_type re;
  fft_type im;
} fft_data;

struct fft_set{
  int N;
  int sgn;
  int factors[64];
  int lf;
  int lt;
  fft_data twiddle[1];
};

fft_object fft_init(int N, int sgn);
void fft_exec(fft_object obj, fft_data *inp, fft_data *oup);
void free_fft(fft_object object);

// CWT
#ifndef cplx_type
#define cplx_type float
#endif

#ifndef cwt_type
#define cwt_type float
#endif

typedef struct cplx_t {
	cplx_type re;
	cplx_type im;
} cplx_data;

typedef struct cwt_set* cwt_object;

struct cwt_set{
  char wave[10];// Wavelet - morl/morlet,paul,dog/dgauss
  int siglength;// Length of Input Data
  int J;// Total Number of Scales
  cwt_type s0;// Smallest scale. It depends on the sampling rate. s0 <= 2 * dt for most wavelets
  cwt_type dt;// Sampling Rate
  cwt_type dj;// Separation between scales. eg., scale = s0 * 2 ^ ( [0:N-1] *dj ) or scale = s0 *[0:N-1] * dj
  char type[10];// Scale Type - Power or Linear
  int pow;// Base of Power in case type = pow. Typical value is pow = 2
  int sflag;
  int pflag;
  int npad;
  int mother;
  cwt_type m;// Wavelet parameter param
  cwt_type smean;// Input Signal mean
  cplx_data *output;
  cwt_type *scale;
  cwt_type *period;
  cwt_type *coi;
  cwt_type params[0];
};

cwt_object cwt_init(const char* wave, cwt_type param, int siglength, cwt_type dt, int J);
void setCWTScales(cwt_object wt, cwt_type s0, cwt_type dj, const char *type, int power);
//void setCWTScaleVector(cwt_object wt, const double *scale, int J, double s0, double dj);
//void setCWTPadding(cwt_object wt, int pad);
void cwavelet(const cwt_type *y, int N, cwt_type dt, int mother, cwt_type param, cwt_type s0, cwt_type dj, int jtot, int npad, cwt_type *wave, cwt_type *scale, cwt_type *period, cwt_type *coi);
void cwt(cwt_object wt, const cwt_type *inp);
void cwt_free(cwt_object object);
void cwt_summary(cwt_object wt);
//void icwt(cwt_object wt, double *cwtop);

//int getCWTScaleLength(int N);
} // namespace sigmap

#endif // CWT_H_
