#include "cwt.h"

namespace sigmap {
// FFT
int dividebyN(int N) {
  while (N % 53 == 0) {
    N = N / 53;
  }
  while (N % 47 == 0) {
    N = N / 47;
  }
  while (N % 43 == 0) {
    N = N / 43;
  }
  while (N % 41 == 0) {
    N = N / 41;
  }
  while (N % 37 == 0) {
    N = N / 37;
  }
  while (N % 31 == 0) {
    N = N / 31;
  }
  while (N % 29 == 0) {
    N = N / 29;
  }
  while (N % 23 == 0) {
    N = N / 23;
  }
  while (N % 17 == 0) {
    N = N / 17;
  }
  while (N % 13 == 0) {
    N = N / 13;
  }
  while (N % 11 == 0) {
    N = N / 11;
  }
  while (N % 8 == 0) {
    N = N / 8;
  }
  while (N % 7 == 0) {
    N = N / 7;
  }
  while (N % 5 == 0) {
    N = N / 5;
  }
  while (N % 4 == 0) {
    N = N / 4;
  }
  while (N % 3 == 0) {
    N = N / 3;
  }
  while (N % 2 == 0) {
    N = N / 2;
  }
  if (N == 1) {
    return 1;
  }
  return 0;
}

int factors(int M, int *arr) {
  int i, N, num, mult, m1, m2;
  i = 0;
  N = M;
  while (N % 53 == 0) {
    N = N / 53;
    arr[i] = 53;
    i++;
  }
  while (N % 47 == 0) {
    N = N / 47;
    arr[i] = 47;
    i++;
  }
  while (N % 43 == 0) {
    N = N / 43;
    arr[i] = 43;
    i++;
  }
  while (N % 41 == 0) {
    N = N / 41;
    arr[i] = 41;
    i++;
  }
  while (N % 37 == 0) {
    N = N / 37;
    arr[i] = 37;
    i++;
  }
  while (N % 31 == 0) {
    N = N / 31;
    arr[i] = 31;
    i++;
  }
  while (N % 29 == 0) {
    N = N / 29;
    arr[i] = 29;
    i++;
  }
  while (N % 23 == 0) {
    N = N / 23;
    arr[i] = 23;
    i++;
  }
  while (N % 19 == 0) {
    N = N / 19;
    arr[i] = 19;
    i++;
  }
  while (N % 17 == 0) {
    N = N / 17;
    arr[i] = 17;
    i++;
  }
  while (N % 13 == 0) {
    N = N / 13;
    arr[i] = 13;
    i++;
  }
  while (N % 11 == 0) {
    N = N / 11;
    arr[i] = 11;
    i++;
  }
  while (N % 8 == 0) {
    N = N / 8;
    arr[i] = 8;
    i++;
  }
  while (N % 7 == 0) {
    N = N / 7;
    arr[i] = 7;
    i++;
  }
  while (N % 5 == 0) {
    N = N / 5;
    arr[i] = 5;
    i++;
  }
  while (N % 4 == 0) {
    N = N / 4;
    arr[i] = 4;
    i++;
  }
  while (N % 3 == 0) {
    N = N / 3;
    arr[i] = 3;
    i++;
  }
  while (N % 2 == 0) {
    N = N / 2;
    arr[i] = 2;
    i++;
  }
  if (N > 31) {
    num = 2;
    while (N > 1) {
      mult = num * 6;
      m1 = mult - 1;
      m2 = mult + 1;
      while (N % m1 == 0) {
        arr[i] = m1;
        i++;
        N = N / m1;
      }
      while (N % m2 == 0) {
        arr[i] = m2;
        i++;
        N = N / m2;
      }
      num += 1;
    }
  }
  return i;
}

void longvectorN(fft_data *sig, int *array, int tx) {
  int L, i, Ls, ct, j, k;
  fft_type theta;
  L = 1;
  ct = 0;
  for (i = 0; i < tx; i++) {
    L = L * array[tx - 1 - i];
    Ls = L / array[tx - 1 - i];
    theta = -1.0 * PI2 / L;
    for (j = 0; j < Ls; j++) {
      for (k = 0; k < array[tx - 1 - i] - 1; k++) {
        sig[ct].re = cos((k + 1) * j * theta);
        sig[ct].im = sin((k + 1) * j * theta);
        ct++;
      }
    }
  }
}

fft_object fft_init(int N, int sgn) {
  fft_object obj = NULL;
  // Change N/2 to N-1 for longvector case
  int twi_len, ct, out;
  out = dividebyN(N);
  if (out == 1) {
    obj =
        (fft_object)malloc(sizeof(struct fft_set) + sizeof(fft_data) * (N - 1));
    obj->lf = factors(N, obj->factors);
    longvectorN(obj->twiddle, obj->factors, obj->lf);
    twi_len = N;
    obj->lt = 0;
  } else {
    int K, M;
    K = (int)pow(2.0, ceil(log10(N) / log10(2.0)));
    if (K < 2 * N - 2) {
      M = K * 2;
    } else {
      M = K;
    }
    obj =
        (fft_object)malloc(sizeof(struct fft_set) + sizeof(fft_data) * (M - 1));
    obj->lf = factors(M, obj->factors);
    longvectorN(obj->twiddle, obj->factors, obj->lf);
    obj->lt = 1;
    twi_len = M;
  }
  obj->N = N;
  obj->sgn = sgn;
  if (sgn == -1) {
    for (ct = 0; ct < twi_len; ct++) {
      (obj->twiddle + ct)->im = -(obj->twiddle + ct)->im;
    }
  }
  return obj;
}

static void bluestein_exp(fft_data *hl, fft_data *hlt, int len, int M) {
  fft_type PI, theta, angle;
  int l2, len2, i;
  PI = 3.1415926535897932384626433832795;
  theta = PI / len;
  l2 = 0;
  len2 = 2 * len;
  for (i = 0; i < len; ++i) {
    angle = theta * l2;
    hlt[i].re = cos(angle);
    hlt[i].im = sin(angle);
    hl[i].re = hlt[i].re;
    hl[i].im = hlt[i].im;
    l2 += 2 * i + 1;
    while (l2 > len2) {
      l2 -= len2;
    }
  }
  for (i = len; i < M - len + 1; i++) {
    hl[i].re = 0.0;
    hl[i].im = 0.0;
  }
  for (i = M - len + 1; i < M; i++) {
    hl[i].re = hlt[M - i].re;
    hl[i].im = hlt[M - i].im;
  }
}

static void bluestein_fft(fft_data *data, fft_data *oup, fft_object obj,
                          int sgn, int N) {
  int K, M, ii, i;
  int def_lt, def_N, def_sgn;
  fft_type scale, temp;
  fft_data *yn;
  fft_data *hk;
  fft_data *tempop;
  fft_data *yno;
  fft_data *hlt;
  obj->lt = 0;
  K = (int)pow(2.0, ceil((double)log10((double)N) / log10((double)2.0)));
  def_lt = 1;
  def_sgn = obj->sgn;
  def_N = obj->N;
  if (K < 2 * N - 2) {
    M = K * 2;
  } else {
    M = K;
  }
  obj->N = M;
  yn = (fft_data *)malloc(sizeof(fft_data) * M);
  hk = (fft_data *)malloc(sizeof(fft_data) * M);
  tempop = (fft_data *)malloc(sizeof(fft_data) * M);
  yno = (fft_data *)malloc(sizeof(fft_data) * M);
  hlt = (fft_data *)malloc(sizeof(fft_data) * N);
  // fft_data* twi = (fft_data*) malloc (sizeof(fft_data) * M);
  bluestein_exp(tempop, hlt, N, M);
  scale = 1.0 / M;
  for (ii = 0; ii < M; ++ii) {
    tempop[ii].im *= scale;
    tempop[ii].re *= scale;
  }
  // fft_object obj = initialize_fft2(M,1);
  fft_exec(obj, tempop, hk);
  if (sgn == 1) {
    for (i = 0; i < N; i++) {
      tempop[i].re = data[i].re * hlt[i].re + data[i].im * hlt[i].im;
      tempop[i].im = -data[i].re * hlt[i].im + data[i].im * hlt[i].re;
    }
  } else {
    for (i = 0; i < N; i++) {
      tempop[i].re = data[i].re * hlt[i].re - data[i].im * hlt[i].im;
      tempop[i].im = data[i].re * hlt[i].im + data[i].im * hlt[i].re;
    }
  }
  for (i = N; i < M; i++) {
    tempop[i].re = 0.0;
    tempop[i].im = 0.0;
  }
  fft_exec(obj, tempop, yn);
  if (sgn == 1) {
    for (i = 0; i < M; i++) {
      temp = yn[i].re * hk[i].re - yn[i].im * hk[i].im;
      yn[i].im = yn[i].re * hk[i].im + yn[i].im * hk[i].re;
      yn[i].re = temp;
    }
  } else {
    for (i = 0; i < M; i++) {
      temp = yn[i].re * hk[i].re + yn[i].im * hk[i].im;
      yn[i].im = -yn[i].re * hk[i].im + yn[i].im * hk[i].re;
      yn[i].re = temp;
    }
  }
  // IFFT
  for (ii = 0; ii < M; ++ii) {
    (obj->twiddle + ii)->im = -(obj->twiddle + ii)->im;
  }
  obj->sgn = -1 * sgn;
  fft_exec(obj, yn, yno);
  if (sgn == 1) {
    for (i = 0; i < N; i++) {
      oup[i].re = yno[i].re * hlt[i].re + yno[i].im * hlt[i].im;
      oup[i].im = -yno[i].re * hlt[i].im + yno[i].im * hlt[i].re;
    }
  } else {
    for (i = 0; i < N; i++) {
      oup[i].re = yno[i].re * hlt[i].re - yno[i].im * hlt[i].im;
      oup[i].im = yno[i].re * hlt[i].im + yno[i].im * hlt[i].re;
    }
  }
  obj->sgn = def_sgn;
  obj->N = def_N;
  obj->lt = def_lt;
  for (ii = 0; ii < M; ++ii) {
    (obj->twiddle + ii)->im = -(obj->twiddle + ii)->im;
  }
  free(yn);
  free(yno);
  free(tempop);
  free(hk);
  free(hlt);
}

static void mixed_radix_dit_rec(fft_data *op, fft_data *ip,
                                const fft_object obj, int sgn, int N, int l,
                                int inc) {
  int radix = 9, m, ll;
  if (N > 1) {
    radix = obj->factors[inc];
    // printf("%d \n",radix);
  }
  if (N == 1) {
    op[0].re = ip[0].re;
    op[0].im = ip[0].im;
  } else if (N == 2) {
    fft_type tau1r, tau1i;
    op[0].re = ip[0].re;
    op[0].im = ip[0].im;
    op[1].re = ip[l].re;
    op[1].im = ip[l].im;
    tau1r = op[0].re;
    tau1i = op[0].im;
    op[0].re = tau1r + op[1].re;
    op[0].im = tau1i + op[1].im;
    op[1].re = tau1r - op[1].re;
    op[1].im = tau1i - op[1].im;
  } else if (N == 3) {
    fft_type tau0r, tau0i, tau1r, tau1i, tau2r, tau2i;
    op[0].re = ip[0].re;
    op[0].im = ip[0].im;
    op[1].re = ip[l].re;
    op[1].im = ip[l].im;
    op[2].re = ip[2 * l].re;
    op[2].im = ip[2 * l].im;
    tau0r = op[1].re + op[2].re;
    tau0i = op[1].im + op[2].im;
    tau1r = sgn * 0.86602540378 * (op[1].re - op[2].re);
    tau1i = sgn * 0.86602540378 * (op[1].im - op[2].im);
    tau2r = op[0].re - tau0r * 0.5000000000;
    tau2i = op[0].im - tau0i * 0.5000000000;
    op[0].re = tau0r + op[0].re;
    op[0].im = tau0i + op[0].im;
    op[1].re = tau2r + tau1i;
    op[1].im = tau2i - tau1r;
    op[2].re = tau2r - tau1i;
    op[2].im = tau2i + tau1r;
    return;
  } else if (N == 4) {
    fft_type tau0r, tau0i, tau1r, tau1i, tau2r, tau2i, tau3r, tau3i;
    op[0].re = ip[0].re;
    op[0].im = ip[0].im;
    op[1].re = ip[l].re;
    op[1].im = ip[l].im;
    op[2].re = ip[2 * l].re;
    op[2].im = ip[2 * l].im;
    op[3].re = ip[3 * l].re;
    op[3].im = ip[3 * l].im;
    tau0r = op[0].re + op[2].re;
    tau0i = op[0].im + op[2].im;
    tau1r = op[0].re - op[2].re;
    tau1i = op[0].im - op[2].im;
    tau2r = op[1].re + op[3].re;
    tau2i = op[1].im + op[3].im;
    tau3r = sgn * (op[1].re - op[3].re);
    tau3i = sgn * (op[1].im - op[3].im);
    op[0].re = tau0r + tau2r;
    op[0].im = tau0i + tau2i;
    op[1].re = tau1r + tau3i;
    op[1].im = tau1i - tau3r;
    op[2].re = tau0r - tau2r;
    op[2].im = tau0i - tau2i;
    op[3].re = tau1r - tau3i;
    op[3].im = tau1i + tau3r;
  } else if (N == 5) {
    fft_type tau0r, tau0i, tau1r, tau1i, tau2r, tau2i, tau3r, tau3i, tau4r,
        tau4i, tau5r, tau5i, tau6r, tau6i;
    fft_type c1, c2, s1, s2;
    op[0].re = ip[0].re;
    op[0].im = ip[0].im;
    op[1].re = ip[l].re;
    op[1].im = ip[l].im;
    op[2].re = ip[2 * l].re;
    op[2].im = ip[2 * l].im;
    op[3].re = ip[3 * l].re;
    op[3].im = ip[3 * l].im;
    op[4].re = ip[4 * l].re;
    op[4].im = ip[4 * l].im;
    c1 = 0.30901699437;
    c2 = -0.80901699437;
    s1 = 0.95105651629;
    s2 = 0.58778525229;
    tau0r = op[1].re + op[4].re;
    tau2r = op[1].re - op[4].re;
    tau0i = op[1].im + op[4].im;
    tau2i = op[1].im - op[4].im;
    tau1r = op[2].re + op[3].re;
    tau3r = op[2].re - op[3].re;
    tau1i = op[2].im + op[3].im;
    tau3i = op[2].im - op[3].im;
    tau4r = c1 * tau0r + c2 * tau1r;
    tau4i = c1 * tau0i + c2 * tau1i;
    // tau5r = sgn * ( s1 * tau2r + s2 * tau3r);
    // tau5i = sgn * ( s1 * tau2i + s2 * tau3i);
    if (sgn == 1) {
      tau5r = s1 * tau2r + s2 * tau3r;
      tau5i = s1 * tau2i + s2 * tau3i;
    } else {
      tau5r = -s1 * tau2r - s2 * tau3r;
      tau5i = -s1 * tau2i - s2 * tau3i;
    }
    tau6r = op[0].re + tau4r;
    tau6i = op[0].im + tau4i;
    op[1].re = tau6r + tau5i;
    op[1].im = tau6i - tau5r;
    op[4].re = tau6r - tau5i;
    op[4].im = tau6i + tau5r;
    tau4r = c2 * tau0r + c1 * tau1r;
    tau4i = c2 * tau0i + c1 * tau1i;
    // tau5r = sgn * ( s2 * tau2r - s1 * tau3r);
    // tau5i = sgn * ( s2 * tau2i - s1 * tau3i);
    if (sgn == 1) {
      tau5r = s2 * tau2r - s1 * tau3r;
      tau5i = s2 * tau2i - s1 * tau3i;
    } else {
      tau5r = -s2 * tau2r + s1 * tau3r;
      tau5i = -s2 * tau2i + s1 * tau3i;
    }
    tau6r = op[0].re + tau4r;
    tau6i = op[0].im + tau4i;
    op[2].re = tau6r + tau5i;
    op[2].im = tau6i - tau5r;
    op[3].re = tau6r - tau5i;
    op[3].im = tau6i + tau5r;
    op[0].re += tau0r + tau1r;
    op[0].im += tau0i + tau1i;
  } else if (N == 7) {
    fft_type tau0r, tau0i, tau1r, tau1i, tau2r, tau2i, tau3r, tau3i, tau4r,
        tau4i, tau5r, tau5i, tau6r, tau6i, tau7r, tau7i;
    fft_type c1, c2, c3, s1, s2, s3;
    op[0].re = ip[0].re;
    op[0].im = ip[0].im;
    op[1].re = ip[l].re;
    op[1].im = ip[l].im;
    op[2].re = ip[2 * l].re;
    op[2].im = ip[2 * l].im;
    op[3].re = ip[3 * l].re;
    op[3].im = ip[3 * l].im;
    op[4].re = ip[4 * l].re;
    op[4].im = ip[4 * l].im;
    op[5].re = ip[5 * l].re;
    op[5].im = ip[5 * l].im;
    op[6].re = ip[6 * l].re;
    op[6].im = ip[6 * l].im;
    c1 = 0.62348980185;
    c2 = -0.22252093395;
    c3 = -0.9009688679;
    s1 = 0.78183148246;
    s2 = 0.97492791218;
    s3 = 0.43388373911;
    tau0r = op[1].re + op[6].re;
    tau3r = op[1].re - op[6].re;
    tau0i = op[1].im + op[6].im;
    tau3i = op[1].im - op[6].im;
    tau1r = op[2].re + op[5].re;
    tau4r = op[2].re - op[5].re;
    tau1i = op[2].im + op[5].im;
    tau4i = op[2].im - op[5].im;
    tau2r = op[3].re + op[4].re;
    tau5r = op[3].re - op[4].re;
    tau2i = op[3].im + op[4].im;
    tau5i = op[3].im - op[4].im;
    tau6r = op[0].re + c1 * tau0r + c2 * tau1r + c3 * tau2r;
    tau6i = op[0].im + c1 * tau0i + c2 * tau1i + c3 * tau2i;
    // tau7r = sgn * ( -s1 * tau3r - s2 * tau4r - s3 * tau5r);
    // tau7i = sgn * ( -s1 * tau3i - s2 * tau4i - s3 * tau5i);
    if (sgn == 1) {
      tau7r = -s1 * tau3r - s2 * tau4r - s3 * tau5r;
      tau7i = -s1 * tau3i - s2 * tau4i - s3 * tau5i;
    } else {
      tau7r = s1 * tau3r + s2 * tau4r + s3 * tau5r;
      tau7i = s1 * tau3i + s2 * tau4i + s3 * tau5i;
    }
    op[1].re = tau6r - tau7i;
    op[6].re = tau6r + tau7i;
    op[1].im = tau6i + tau7r;
    op[6].im = tau6i - tau7r;
    tau6r = op[0].re + c2 * tau0r + c3 * tau1r + c1 * tau2r;
    tau6i = op[0].im + c2 * tau0i + c3 * tau1i + c1 * tau2i;
    // tau7r = sgn * ( -s2 * tau3r + s3 * tau4r + s1 * tau5r);
    // tau7i = sgn * ( -s2 * tau3i + s3 * tau4i + s1 * tau5i);
    if (sgn == 1) {
      tau7r = -s2 * tau3r + s3 * tau4r + s1 * tau5r;
      tau7i = -s2 * tau3i + s3 * tau4i + s1 * tau5i;
    } else {
      tau7r = s2 * tau3r - s3 * tau4r - s1 * tau5r;
      tau7i = s2 * tau3i - s3 * tau4i - s1 * tau5i;
    }
    op[2].re = tau6r - tau7i;
    op[5].re = tau6r + tau7i;
    op[2].im = tau6i + tau7r;
    op[5].im = tau6i - tau7r;
    tau6r = op[0].re + c3 * tau0r + c1 * tau1r + c2 * tau2r;
    tau6i = op[0].im + c3 * tau0i + c1 * tau1i + c2 * tau2i;
    // tau7r = sgn * ( -s3 * tau3r + s1 * tau4r - s2 * tau5r);
    // tau7i = sgn * ( -s3 * tau3i + s1 * tau4i - s2 * tau5i);
    if (sgn == 1) {
      tau7r = -s3 * tau3r + s1 * tau4r - s2 * tau5r;
      tau7i = -s3 * tau3i + s1 * tau4i - s2 * tau5i;
    } else {
      tau7r = s3 * tau3r - s1 * tau4r + s2 * tau5r;
      tau7i = s3 * tau3i - s1 * tau4i + s2 * tau5i;
    }
    op[3].re = tau6r - tau7i;
    op[4].re = tau6r + tau7i;
    op[3].im = tau6i + tau7r;
    op[4].im = tau6i - tau7r;
    op[0].re += tau0r + tau1r + tau2r;
    op[0].im += tau0i + tau1i + tau2i;
  } else if (N == 8) {
    fft_type tau0r, tau0i, tau1r, tau1i, tau2r, tau2i, tau3r, tau3i, tau4r,
        tau4i, tau5r, tau5i, tau6r, tau6i, tau7r, tau7i, tau8r, tau8i, tau9r,
        tau9i;
    fft_type c1, s1, temp1r, temp1i, temp2r, temp2i;
    op[0].re = ip[0].re;
    op[0].im = ip[0].im;
    op[1].re = ip[l].re;
    op[1].im = ip[l].im;
    op[2].re = ip[2 * l].re;
    op[2].im = ip[2 * l].im;
    op[3].re = ip[3 * l].re;
    op[3].im = ip[3 * l].im;
    op[4].re = ip[4 * l].re;
    op[4].im = ip[4 * l].im;
    op[5].re = ip[5 * l].re;
    op[5].im = ip[5 * l].im;
    op[6].re = ip[6 * l].re;
    op[6].im = ip[6 * l].im;
    op[7].re = ip[7 * l].re;
    op[7].im = ip[7 * l].im;
    c1 = 0.70710678118654752440084436210485;
    s1 = 0.70710678118654752440084436210485;
    tau0r = op[0].re + op[4].re;
    tau4r = op[0].re - op[4].re;
    tau0i = op[0].im + op[4].im;
    tau4i = op[0].im - op[4].im;
    tau1r = op[1].re + op[7].re;
    tau5r = op[1].re - op[7].re;
    tau1i = op[1].im + op[7].im;
    tau5i = op[1].im - op[7].im;
    tau2r = op[3].re + op[5].re;
    tau6r = op[3].re - op[5].re;
    tau2i = op[3].im + op[5].im;
    tau6i = op[3].im - op[5].im;
    tau3r = op[2].re + op[6].re;
    tau7r = op[2].re - op[6].re;
    tau3i = op[2].im + op[6].im;
    tau7i = op[2].im - op[6].im;
    op[0].re = tau0r + tau1r + tau2r + tau3r;
    op[0].im = tau0i + tau1i + tau2i + tau3i;
    op[4].re = tau0r - tau1r - tau2r + tau3r;
    op[4].im = tau0i - tau1i - tau2i + tau3i;
    temp1r = tau1r - tau2r;
    temp1i = tau1i - tau2i;
    temp2r = tau5r + tau6r;
    temp2i = tau5i + tau6i;
    tau8r = tau4r + c1 * temp1r;
    tau8i = tau4i + c1 * temp1i;
    // tau9r = sgn * ( -s1 * temp2r - tau7r);
    // tau9i = sgn * ( -s1 * temp2i - tau7i);
    if (sgn == 1) {
      tau9r = -s1 * temp2r - tau7r;
      tau9i = -s1 * temp2i - tau7i;
    } else {
      tau9r = s1 * temp2r + tau7r;
      tau9i = s1 * temp2i + tau7i;
    }
    op[1].re = tau8r - tau9i;
    op[1].im = tau8i + tau9r;
    op[7].re = tau8r + tau9i;
    op[7].im = tau8i - tau9r;
    tau8r = tau0r - tau3r;
    tau8i = tau0i - tau3i;
    // tau9r = sgn * ( -tau5r + tau6r);
    // tau9i = sgn * ( -tau5i + tau6i);
    if (sgn == 1) {
      tau9r = -tau5r + tau6r;
      tau9i = -tau5i + tau6i;
    } else {
      tau9r = tau5r - tau6r;
      tau9i = tau5i - tau6i;
    }
    op[2].re = tau8r - tau9i;
    op[2].im = tau8i + tau9r;
    op[6].re = tau8r + tau9i;
    op[6].im = tau8i - tau9r;
    tau8r = tau4r - c1 * temp1r;
    tau8i = tau4i - c1 * temp1i;
    // tau9r = sgn * ( -s1 * temp2r + tau7r);
    // tau9i = sgn * ( -s1 * temp2i + tau7i);
    if (sgn == 1) {
      tau9r = -s1 * temp2r + tau7r;
      tau9i = -s1 * temp2i + tau7i;
    } else {
      tau9r = s1 * temp2r - tau7r;
      tau9i = s1 * temp2i - tau7i;
    }
    op[3].re = tau8r - tau9i;
    op[3].im = tau8i + tau9r;
    op[5].re = tau8r + tau9i;
    op[5].im = tau8i - tau9r;
  } else if (radix == 2) {
    int k, tkm1, ind;
    fft_type wlr, wli;
    fft_type tau1r, tau1i, tau2r, tau2i;
    m = N / 2;
    ll = 2 * l;
    mixed_radix_dit_rec(op, ip, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + m, ip + l, obj, sgn, m, ll, inc + 1);
    for (k = 0; k < m; k++) {
      ind = m - 1 + k;
      wlr = (obj->twiddle + ind)->re;
      wli = (obj->twiddle + ind)->im;
      tkm1 = k + m;
      tau1r = op[k].re;
      tau1i = op[k].im;
      tau2r = op[tkm1].re * wlr - op[tkm1].im * wli;
      tau2i = op[tkm1].im * wlr + op[tkm1].re * wli;
      op[k].re = tau1r + tau2r;
      op[k].im = tau1i + tau2i;
      op[tkm1].re = tau1r - tau2r;
      op[tkm1].im = tau1i - tau2i;
    }
  } else if (radix == 3) {
    int k, tkm1, tkm2, ind;
    fft_type wlr, wli, wl2r, wl2i;
    fft_type tau0r, tau0i, tau1r, tau1i, tau2r, tau2i;
    fft_type ar, ai, br, bi, cr, ci;
    m = N / 3;
    ll = 3 * l;
    mixed_radix_dit_rec(op, ip, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + m, ip + l, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + 2 * m, ip + 2 * l, obj, sgn, m, ll, inc + 1);
    // printf("%d \n",inc);
    // mixed_radix3_dit_rec(op,ip,obj,sgn,ll,m);
    for (k = 0; k < m; ++k) {
      ind = m - 1 + 2 * k;
      wlr = (obj->twiddle + ind)->re;
      wli = (obj->twiddle + ind)->im;
      ind++;
      wl2r = (obj->twiddle + ind)->re;
      wl2i = (obj->twiddle + ind)->im;
      tkm1 = k + m;
      tkm2 = tkm1 + m;
      ar = op[k].re;
      ai = op[k].im;
      br = op[tkm1].re * wlr - op[tkm1].im * wli;
      bi = op[tkm1].im * wlr + op[tkm1].re * wli;
      cr = op[tkm2].re * wl2r - op[tkm2].im * wl2i;
      ci = op[tkm2].im * wl2r + op[tkm2].re * wl2i;
      tau0r = br + cr;
      tau0i = bi + ci;
      tau1r = sgn * 0.86602540378 * (br - cr);
      tau1i = sgn * 0.86602540378 * (bi - ci);
      tau2r = ar - tau0r * 0.5000000000;
      tau2i = ai - tau0i * 0.5000000000;
      op[k].re = ar + tau0r;
      op[k].im = ai + tau0i;
      op[tkm1].re = tau2r + tau1i;
      op[tkm1].im = tau2i - tau1r;
      op[tkm2].re = tau2r - tau1i;
      op[tkm2].im = tau2i + tau1r;
    }
  } else if (radix == 4) {
    int k, tkm1, tkm2, tkm3, ind;
    fft_type wlr, wli, wl2r, wl2i, wl3r, wl3i;
    fft_type tau0r, tau0i, tau1r, tau1i, tau2r, tau2i, tau3r, tau3i;
    fft_type ar, ai, br, bi, cr, ci, dr, di;
    m = N / 4;
    ll = 4 * l;
    mixed_radix_dit_rec(op, ip, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + m, ip + l, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + 2 * m, ip + 2 * l, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + 3 * m, ip + 3 * l, obj, sgn, m, ll, inc + 1);
    // mixed_radix4_dit_rec(op,ip,obj,sgn,ll,m);
    tkm1 = m;
    tkm2 = tkm1 + m;
    tkm3 = tkm2 + m;
    ar = op[0].re;
    ai = op[0].im;
    br = op[tkm1].re;
    bi = op[tkm1].im;
    cr = op[tkm2].re;
    ci = op[tkm2].im;
    dr = op[tkm3].re;
    di = op[tkm3].im;
    tau0r = ar + cr;
    tau0i = ai + ci;
    tau1r = ar - cr;
    tau1i = ai - ci;
    tau2r = br + dr;
    tau2i = bi + di;
    tau3r = sgn * (br - dr);
    tau3i = sgn * (bi - di);
    op[0].re = tau0r + tau2r;
    op[0].im = tau0i + tau2i;
    op[tkm1].re = tau1r + tau3i;
    op[tkm1].im = tau1i - tau3r;
    op[tkm2].re = tau0r - tau2r;
    op[tkm2].im = tau0i - tau2i;
    op[tkm3].re = tau1r - tau3i;
    op[tkm3].im = tau1i + tau3r;
    for (k = 1; k < m; k++) {
      ind = m - 1 + 3 * k;
      wlr = (obj->twiddle + ind)->re;
      wli = (obj->twiddle + ind)->im;
      ind++;
      wl2r = (obj->twiddle + ind)->re;
      wl2i = (obj->twiddle + ind)->im;
      ind++;
      wl3r = (obj->twiddle + ind)->re;
      wl3i = (obj->twiddle + ind)->im;
      tkm1 = k + m;
      tkm2 = tkm1 + m;
      tkm3 = tkm2 + m;
      ar = op[k].re;
      ai = op[k].im;
      br = op[tkm1].re * wlr - op[tkm1].im * wli;
      bi = op[tkm1].im * wlr + op[tkm1].re * wli;
      cr = op[tkm2].re * wl2r - op[tkm2].im * wl2i;
      ci = op[tkm2].im * wl2r + op[tkm2].re * wl2i;
      dr = op[tkm3].re * wl3r - op[tkm3].im * wl3i;
      di = op[tkm3].im * wl3r + op[tkm3].re * wl3i;
      tau0r = ar + cr;
      tau0i = ai + ci;
      tau1r = ar - cr;
      tau1i = ai - ci;
      tau2r = br + dr;
      tau2i = bi + di;
      tau3r = sgn * (br - dr);
      tau3i = sgn * (bi - di);
      op[k].re = tau0r + tau2r;
      op[k].im = tau0i + tau2i;
      op[tkm1].re = tau1r + tau3i;
      op[tkm1].im = tau1i - tau3r;
      op[tkm2].re = tau0r - tau2r;
      op[tkm2].im = tau0i - tau2i;
      op[tkm3].re = tau1r - tau3i;
      op[tkm3].im = tau1i + tau3r;
    }
  } else if (radix == 5) {
    int k, tkm1, tkm2, tkm3, tkm4, ind;
    fft_type wlr, wli, wl2r, wl2i, wl3r, wl3i, wl4r, wl4i;
    fft_type tau0r, tau0i, tau1r, tau1i, tau2r, tau2i, tau3r, tau3i;
    fft_type ar, ai, br, bi, cr, ci, dr, di, er, ei;
    fft_type tau4r, tau4i, tau5r, tau5i, tau6r, tau6i;
    fft_type c1, c2, s1, s2;
    m = N / 5;
    ll = 5 * l;
    mixed_radix_dit_rec(op, ip, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + m, ip + l, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + 2 * m, ip + 2 * l, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + 3 * m, ip + 3 * l, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + 4 * m, ip + 4 * l, obj, sgn, m, ll, inc + 1);
    // printf("%d \n",inc);
    // mixed_radix3_dit_rec(op,ip,obj,sgn,ll,m);
    c1 = 0.30901699437;
    c2 = -0.80901699437;
    s1 = 0.95105651629;
    s2 = 0.58778525229;
    tkm1 = m;
    tkm2 = tkm1 + m;
    tkm3 = tkm2 + m;
    tkm4 = tkm3 + m;
    ar = op[0].re;
    ai = op[0].im;
    br = op[tkm1].re;
    bi = op[tkm1].im;
    cr = op[tkm2].re;
    ci = op[tkm2].im;
    dr = op[tkm3].re;
    di = op[tkm3].im;
    er = op[tkm4].re;
    ei = op[tkm4].im;
    tau0r = br + er;
    tau0i = bi + ei;
    tau1r = cr + dr;
    tau1i = ci + di;
    tau2r = br - er;
    tau2i = bi - ei;
    tau3r = cr - dr;
    tau3i = ci - di;
    op[0].re = ar + tau0r + tau1r;
    op[0].im = ai + tau0i + tau1i;
    tau4r = c1 * tau0r + c2 * tau1r;
    tau4i = c1 * tau0i + c2 * tau1i;
    tau5r = sgn * (s1 * tau2r + s2 * tau3r);
    tau5i = sgn * (s1 * tau2i + s2 * tau3i);
    tau6r = ar + tau4r;
    tau6i = ai + tau4i;
    op[tkm1].re = tau6r + tau5i;
    op[tkm1].im = tau6i - tau5r;
    op[tkm4].re = tau6r - tau5i;
    op[tkm4].im = tau6i + tau5r;
    tau4r = c2 * tau0r + c1 * tau1r;
    tau4i = c2 * tau0i + c1 * tau1i;
    tau5r = sgn * (s2 * tau2r - s1 * tau3r);
    tau5i = sgn * (s2 * tau2i - s1 * tau3i);
    tau6r = ar + tau4r;
    tau6i = ai + tau4i;
    op[tkm2].re = tau6r + tau5i;
    op[tkm2].im = tau6i - tau5r;
    op[tkm3].re = tau6r - tau5i;
    op[tkm3].im = tau6i + tau5r;
    for (k = 1; k < m; k++) {
      ind = m - 1 + 4 * k;
      wlr = (obj->twiddle + ind)->re;
      wli = (obj->twiddle + ind)->im;
      ind++;
      wl2r = (obj->twiddle + ind)->re;
      wl2i = (obj->twiddle + ind)->im;
      ind++;
      wl3r = (obj->twiddle + ind)->re;
      wl3i = (obj->twiddle + ind)->im;
      ind++;
      wl4r = (obj->twiddle + ind)->re;
      wl4i = (obj->twiddle + ind)->im;
      tkm1 = k + m;
      tkm2 = tkm1 + m;
      tkm3 = tkm2 + m;
      tkm4 = tkm3 + m;
      ar = op[k].re;
      ai = op[k].im;
      br = op[tkm1].re * wlr - op[tkm1].im * wli;
      bi = op[tkm1].im * wlr + op[tkm1].re * wli;
      cr = op[tkm2].re * wl2r - op[tkm2].im * wl2i;
      ci = op[tkm2].im * wl2r + op[tkm2].re * wl2i;
      dr = op[tkm3].re * wl3r - op[tkm3].im * wl3i;
      di = op[tkm3].im * wl3r + op[tkm3].re * wl3i;
      er = op[tkm4].re * wl4r - op[tkm4].im * wl4i;
      ei = op[tkm4].im * wl4r + op[tkm4].re * wl4i;
      tau0r = br + er;
      tau0i = bi + ei;
      tau1r = cr + dr;
      tau1i = ci + di;
      tau2r = br - er;
      tau2i = bi - ei;
      tau3r = cr - dr;
      tau3i = ci - di;
      op[k].re = ar + tau0r + tau1r;
      op[k].im = ai + tau0i + tau1i;
      tau4r = c1 * tau0r + c2 * tau1r;
      tau4i = c1 * tau0i + c2 * tau1i;
      // tau5r = sgn * ( s1 * tau2r + s2 * tau3r);
      // tau5i = sgn * ( s1 * tau2i + s2 * tau3i);
      if (sgn == 1) {
        tau5r = s1 * tau2r + s2 * tau3r;
        tau5i = s1 * tau2i + s2 * tau3i;
      } else {
        tau5r = -s1 * tau2r - s2 * tau3r;
        tau5i = -s1 * tau2i - s2 * tau3i;
      }
      tau6r = ar + tau4r;
      tau6i = ai + tau4i;
      op[tkm1].re = tau6r + tau5i;
      op[tkm1].im = tau6i - tau5r;
      op[tkm4].re = tau6r - tau5i;
      op[tkm4].im = tau6i + tau5r;
      tau4r = c2 * tau0r + c1 * tau1r;
      tau4i = c2 * tau0i + c1 * tau1i;
      // tau5r = sgn * ( s2 * tau2r - s1 * tau3r);
      // tau5i = sgn * ( s2 * tau2i - s1 * tau3i);
      if (sgn == 1) {
        tau5r = s2 * tau2r - s1 * tau3r;
        tau5i = s2 * tau2i - s1 * tau3i;
      } else {
        tau5r = -s2 * tau2r + s1 * tau3r;
        tau5i = -s2 * tau2i + s1 * tau3i;
      }
      tau6r = ar + tau4r;
      tau6i = ai + tau4i;
      op[tkm2].re = tau6r + tau5i;
      op[tkm2].im = tau6i - tau5r;
      op[tkm3].re = tau6r - tau5i;
      op[tkm3].im = tau6i + tau5r;
    }
  } else if (radix == 7) {
    int k, tkm1, tkm2, tkm3, tkm4, tkm5, tkm6, ind;
    fft_type wlr, wli, wl2r, wl2i, wl3r, wl3i, wl4r, wl4i, wl5r, wl5i, wl6r,
        wl6i;
    fft_type tau0r, tau0i, tau1r, tau1i, tau2r, tau2i, tau3r, tau3i;
    fft_type ar, ai, br, bi, cr, ci, dr, di, er, ei, fr, fi, gr, gi;
    fft_type tau4r, tau4i, tau5r, tau5i, tau6r, tau6i, tau7r, tau7i;
    fft_type c1, c2, c3, s1, s2, s3;
    m = N / 7;
    ll = 7 * l;
    mixed_radix_dit_rec(op, ip, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + m, ip + l, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + 2 * m, ip + 2 * l, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + 3 * m, ip + 3 * l, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + 4 * m, ip + 4 * l, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + 5 * m, ip + 5 * l, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + 6 * m, ip + 6 * l, obj, sgn, m, ll, inc + 1);
    // printf("%d \n",inc);
    // mixed_radix3_dit_rec(op,ip,obj,sgn,ll,m);
    c1 = 0.62348980185;
    c2 = -0.22252093395;
    c3 = -0.9009688679;
    s1 = 0.78183148246;
    s2 = 0.97492791218;
    s3 = 0.43388373911;
    tkm1 = m;
    tkm2 = tkm1 + m;
    tkm3 = tkm2 + m;
    tkm4 = tkm3 + m;
    tkm5 = tkm4 + m;
    tkm6 = tkm5 + m;
    ar = op[0].re;
    ai = op[0].im;
    br = op[tkm1].re;
    bi = op[tkm1].im;
    cr = op[tkm2].re;
    ci = op[tkm2].im;
    dr = op[tkm3].re;
    di = op[tkm3].im;
    er = op[tkm4].re;
    ei = op[tkm4].im;
    fr = op[tkm5].re;
    fi = op[tkm5].im;
    gr = op[tkm6].re;
    gi = op[tkm6].im;
    tau0r = br + gr;
    tau3r = br - gr;
    tau0i = bi + gi;
    tau3i = bi - gi;
    tau1r = cr + fr;
    tau4r = cr - fr;
    tau1i = ci + fi;
    tau4i = ci - fi;
    tau2r = dr + er;
    tau5r = dr - er;
    tau2i = di + ei;
    tau5i = di - ei;
    op[0].re = ar + tau0r + tau1r + tau2r;
    op[0].im = ai + tau0i + tau1i + tau2i;
    tau6r = ar + c1 * tau0r + c2 * tau1r + c3 * tau2r;
    tau6i = ai + c1 * tau0i + c2 * tau1i + c3 * tau2i;
    // tau7r = sgn * ( -s1 * tau3r - s2 * tau4r - s3 * tau5r);
    // tau7i = sgn * ( -s1 * tau3i - s2 * tau4i - s3 * tau5i);
    if (sgn == 1) {
      tau7r = -s1 * tau3r - s2 * tau4r - s3 * tau5r;
      tau7i = -s1 * tau3i - s2 * tau4i - s3 * tau5i;
    } else {
      tau7r = s1 * tau3r + s2 * tau4r + s3 * tau5r;
      tau7i = s1 * tau3i + s2 * tau4i + s3 * tau5i;
    }
    op[tkm1].re = tau6r - tau7i;
    op[tkm1].im = tau6i + tau7r;
    op[tkm6].re = tau6r + tau7i;
    op[tkm6].im = tau6i - tau7r;
    tau6r = ar + c2 * tau0r + c3 * tau1r + c1 * tau2r;
    tau6i = ai + c2 * tau0i + c3 * tau1i + c1 * tau2i;
    // tau7r = sgn * ( -s2 * tau3r + s3 * tau4r + s1 * tau5r);
    // tau7i = sgn * ( -s2 * tau3i + s3 * tau4i + s1 * tau5i);
    if (sgn == 1) {
      tau7r = -s2 * tau3r + s3 * tau4r + s1 * tau5r;
      tau7i = -s2 * tau3i + s3 * tau4i + s1 * tau5i;
    } else {
      tau7r = s2 * tau3r - s3 * tau4r - s1 * tau5r;
      tau7i = s2 * tau3i - s3 * tau4i - s1 * tau5i;
    }
    op[tkm2].re = tau6r - tau7i;
    op[tkm2].im = tau6i + tau7r;
    op[tkm5].re = tau6r + tau7i;
    op[tkm5].im = tau6i - tau7r;
    tau6r = ar + c3 * tau0r + c1 * tau1r + c2 * tau2r;
    tau6i = ai + c3 * tau0i + c1 * tau1i + c2 * tau2i;
    // tau7r = sgn * ( -s3 * tau3r + s1 * tau4r - s2 * tau5r);
    // tau7i = sgn * ( -s3 * tau3i + s1 * tau4i - s2 * tau5i);
    if (sgn == 1) {
      tau7r = -s3 * tau3r + s1 * tau4r - s2 * tau5r;
      tau7i = -s3 * tau3i + s1 * tau4i - s2 * tau5i;
    } else {
      tau7r = s3 * tau3r - s1 * tau4r + s2 * tau5r;
      tau7i = s3 * tau3i - s1 * tau4i + s2 * tau5i;
    }
    op[tkm3].re = tau6r - tau7i;
    op[tkm3].im = tau6i + tau7r;
    op[tkm4].re = tau6r + tau7i;
    op[tkm4].im = tau6i - tau7r;
    for (k = 1; k < m; k++) {
      ind = m - 1 + 6 * k;
      wlr = (obj->twiddle + ind)->re;
      wli = (obj->twiddle + ind)->im;
      ind++;
      wl2r = (obj->twiddle + ind)->re;
      wl2i = (obj->twiddle + ind)->im;
      ind++;
      wl3r = (obj->twiddle + ind)->re;
      wl3i = (obj->twiddle + ind)->im;
      ind++;
      wl4r = (obj->twiddle + ind)->re;
      wl4i = (obj->twiddle + ind)->im;
      ind++;
      wl5r = (obj->twiddle + ind)->re;
      wl5i = (obj->twiddle + ind)->im;
      ind++;
      wl6r = (obj->twiddle + ind)->re;
      wl6i = (obj->twiddle + ind)->im;
      tkm1 = k + m;
      tkm2 = tkm1 + m;
      tkm3 = tkm2 + m;
      tkm4 = tkm3 + m;
      tkm5 = tkm4 + m;
      tkm6 = tkm5 + m;
      ar = op[k].re;
      ai = op[k].im;
      br = op[tkm1].re * wlr - op[tkm1].im * wli;
      bi = op[tkm1].im * wlr + op[tkm1].re * wli;
      cr = op[tkm2].re * wl2r - op[tkm2].im * wl2i;
      ci = op[tkm2].im * wl2r + op[tkm2].re * wl2i;
      dr = op[tkm3].re * wl3r - op[tkm3].im * wl3i;
      di = op[tkm3].im * wl3r + op[tkm3].re * wl3i;
      er = op[tkm4].re * wl4r - op[tkm4].im * wl4i;
      ei = op[tkm4].im * wl4r + op[tkm4].re * wl4i;
      fr = op[tkm5].re * wl5r - op[tkm5].im * wl5i;
      fi = op[tkm5].im * wl5r + op[tkm5].re * wl5i;
      gr = op[tkm6].re * wl6r - op[tkm6].im * wl6i;
      gi = op[tkm6].im * wl6r + op[tkm6].re * wl6i;
      tau0r = br + gr;
      tau3r = br - gr;
      tau0i = bi + gi;
      tau3i = bi - gi;
      tau1r = cr + fr;
      tau4r = cr - fr;
      tau1i = ci + fi;
      tau4i = ci - fi;
      tau2r = dr + er;
      tau5r = dr - er;
      tau2i = di + ei;
      tau5i = di - ei;
      op[k].re = ar + tau0r + tau1r + tau2r;
      op[k].im = ai + tau0i + tau1i + tau2i;
      tau6r = ar + c1 * tau0r + c2 * tau1r + c3 * tau2r;
      tau6i = ai + c1 * tau0i + c2 * tau1i + c3 * tau2i;
      // tau7r = sgn * ( -s1 * tau3r - s2 * tau4r - s3 * tau5r);
      // tau7i = sgn * ( -s1 * tau3i - s2 * tau4i - s3 * tau5i);
      if (sgn == 1) {
        tau7r = -s1 * tau3r - s2 * tau4r - s3 * tau5r;
        tau7i = -s1 * tau3i - s2 * tau4i - s3 * tau5i;
      } else {
        tau7r = s1 * tau3r + s2 * tau4r + s3 * tau5r;
        tau7i = s1 * tau3i + s2 * tau4i + s3 * tau5i;
      }
      op[tkm1].re = tau6r - tau7i;
      op[tkm1].im = tau6i + tau7r;
      op[tkm6].re = tau6r + tau7i;
      op[tkm6].im = tau6i - tau7r;
      tau6r = ar + c2 * tau0r + c3 * tau1r + c1 * tau2r;
      tau6i = ai + c2 * tau0i + c3 * tau1i + c1 * tau2i;
      // tau7r = sgn * ( -s2 * tau3r + s3 * tau4r + s1 * tau5r);
      // tau7i = sgn * ( -s2 * tau3i + s3 * tau4i + s1 * tau5i);
      if (sgn == 1) {
        tau7r = -s2 * tau3r + s3 * tau4r + s1 * tau5r;
        tau7i = -s2 * tau3i + s3 * tau4i + s1 * tau5i;
      } else {
        tau7r = s2 * tau3r - s3 * tau4r - s1 * tau5r;
        tau7i = s2 * tau3i - s3 * tau4i - s1 * tau5i;
      }
      op[tkm2].re = tau6r - tau7i;
      op[tkm2].im = tau6i + tau7r;
      op[tkm5].re = tau6r + tau7i;
      op[tkm5].im = tau6i - tau7r;
      tau6r = ar + c3 * tau0r + c1 * tau1r + c2 * tau2r;
      tau6i = ai + c3 * tau0i + c1 * tau1i + c2 * tau2i;
      // tau7r = sgn * ( -s3 * tau3r + s1 * tau4r - s2 * tau5r);
      // tau7i = sgn * ( -s3 * tau3i + s1 * tau4i - s2 * tau5i);
      if (sgn == 1) {
        tau7r = -s3 * tau3r + s1 * tau4r - s2 * tau5r;
        tau7i = -s3 * tau3i + s1 * tau4i - s2 * tau5i;
      } else {
        tau7r = s3 * tau3r - s1 * tau4r + s2 * tau5r;
        tau7i = s3 * tau3i - s1 * tau4i + s2 * tau5i;
      }
      op[tkm3].re = tau6r - tau7i;
      op[tkm3].im = tau6i + tau7r;
      op[tkm4].re = tau6r + tau7i;
      op[tkm4].im = tau6i - tau7r;
    }
  } else if (radix == 8) {
    int k, tkm1, tkm2, tkm3, tkm4, tkm5, tkm6, tkm7, ind;
    fft_type wlr, wli, wl2r, wl2i, wl3r, wl3i, wl4r, wl4i, wl5r, wl5i, wl6r,
        wl6i, wl7r, wl7i;
    fft_type tau0r, tau0i, tau1r, tau1i, tau2r, tau2i, tau3r, tau3i;
    fft_type ar, ai, br, bi, cr, ci, dr, di, er, ei, fr, fi, gr, gi, hr, hi;
    fft_type tau4r, tau4i, tau5r, tau5i, tau6r, tau6i, tau7r, tau7i, tau8r,
        tau8i, tau9r, tau9i;
    fft_type c1, s1, temp1r, temp1i, temp2r, temp2i;
    m = N / 8;
    ll = 8 * l;
    mixed_radix_dit_rec(op, ip, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + m, ip + l, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + 2 * m, ip + 2 * l, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + 3 * m, ip + 3 * l, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + 4 * m, ip + 4 * l, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + 5 * m, ip + 5 * l, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + 6 * m, ip + 6 * l, obj, sgn, m, ll, inc + 1);
    mixed_radix_dit_rec(op + 7 * m, ip + 7 * l, obj, sgn, m, ll, inc + 1);
    // printf("%d \n",inc);
    // mixed_radix3_dit_rec(op,ip,obj,sgn,ll,m);
    c1 = 0.70710678118654752440084436210485;
    s1 = 0.70710678118654752440084436210485;
    for (k = 0; k < m; k++) {
      ind = m - 1 + 7 * k;
      wlr = (obj->twiddle + ind)->re;
      wli = (obj->twiddle + ind)->im;
      ind++;
      wl2r = (obj->twiddle + ind)->re;
      wl2i = (obj->twiddle + ind)->im;
      ind++;
      wl3r = (obj->twiddle + ind)->re;
      wl3i = (obj->twiddle + ind)->im;
      ind++;
      wl4r = (obj->twiddle + ind)->re;
      wl4i = (obj->twiddle + ind)->im;
      ind++;
      wl5r = (obj->twiddle + ind)->re;
      wl5i = (obj->twiddle + ind)->im;
      ind++;
      wl6r = (obj->twiddle + ind)->re;
      wl6i = (obj->twiddle + ind)->im;
      ind++;
      wl7r = (obj->twiddle + ind)->re;
      wl7i = (obj->twiddle + ind)->im;
      tkm1 = k + m;
      tkm2 = tkm1 + m;
      tkm3 = tkm2 + m;
      tkm4 = tkm3 + m;
      tkm5 = tkm4 + m;
      tkm6 = tkm5 + m;
      tkm7 = tkm6 + m;
      ar = op[k].re;
      ai = op[k].im;
      br = op[tkm1].re * wlr - op[tkm1].im * wli;
      bi = op[tkm1].im * wlr + op[tkm1].re * wli;
      cr = op[tkm2].re * wl2r - op[tkm2].im * wl2i;
      ci = op[tkm2].im * wl2r + op[tkm2].re * wl2i;
      dr = op[tkm3].re * wl3r - op[tkm3].im * wl3i;
      di = op[tkm3].im * wl3r + op[tkm3].re * wl3i;
      er = op[tkm4].re * wl4r - op[tkm4].im * wl4i;
      ei = op[tkm4].im * wl4r + op[tkm4].re * wl4i;
      fr = op[tkm5].re * wl5r - op[tkm5].im * wl5i;
      fi = op[tkm5].im * wl5r + op[tkm5].re * wl5i;
      gr = op[tkm6].re * wl6r - op[tkm6].im * wl6i;
      gi = op[tkm6].im * wl6r + op[tkm6].re * wl6i;
      hr = op[tkm7].re * wl7r - op[tkm7].im * wl7i;
      hi = op[tkm7].im * wl7r + op[tkm7].re * wl7i;
      tau0r = ar + er;
      tau4r = ar - er;
      tau0i = ai + ei;
      tau4i = ai - ei;
      tau1r = br + hr;
      tau5r = br - hr;
      tau1i = bi + hi;
      tau5i = bi - hi;
      tau2r = dr + fr;
      tau6r = dr - fr;
      tau6i = di - fi;
      tau2i = di + fi;
      tau3r = cr + gr;
      tau7r = cr - gr;
      tau7i = ci - gi;
      tau3i = ci + gi;
      op[k].re = tau0r + tau1r + tau2r + tau3r;
      op[k].im = tau0i + tau1i + tau2i + tau3i;
      op[tkm4].re = tau0r - tau1r - tau2r + tau3r;
      op[tkm4].im = tau0i - tau1i - tau2i + tau3i;
      temp1r = tau1r - tau2r;
      temp1i = tau1i - tau2i;
      temp2r = tau5r + tau6r;
      temp2i = tau5i + tau6i;
      tau8r = tau4r + c1 * temp1r;
      tau8i = tau4i + c1 * temp1i;
      // tau9r = sgn * ( -s1 * temp2r - tau7r);
      // tau9i = sgn * ( -s1 * temp2i - tau7i);
      if (sgn == 1) {
        tau9r = -s1 * temp2r - tau7r;
        tau9i = -s1 * temp2i - tau7i;
      } else {
        tau9r = s1 * temp2r + tau7r;
        tau9i = s1 * temp2i + tau7i;
      }
      op[tkm1].re = tau8r - tau9i;
      op[tkm1].im = tau8i + tau9r;
      op[tkm7].re = tau8r + tau9i;
      op[tkm7].im = tau8i - tau9r;
      tau8r = tau0r - tau3r;
      tau8i = tau0i - tau3i;
      // tau9r = sgn * ( -tau5r + tau6r);
      // tau9i = sgn * ( -tau5i + tau6i);
      if (sgn == 1) {
        tau9r = -tau5r + tau6r;
        tau9i = -tau5i + tau6i;
      } else {
        tau9r = tau5r - tau6r;
        tau9i = tau5i - tau6i;
      }
      op[tkm2].re = tau8r - tau9i;
      op[tkm2].im = tau8i + tau9r;
      op[tkm6].re = tau8r + tau9i;
      op[tkm6].im = tau8i - tau9r;
      tau8r = tau4r - c1 * temp1r;
      tau8i = tau4i - c1 * temp1i;
      // tau9r = sgn * ( -s1 * temp2r + tau7r);
      // tau9i = sgn * ( -s1 * temp2i + tau7i);
      if (sgn == 1) {
        tau9r = -s1 * temp2r + tau7r;
        tau9i = -s1 * temp2i + tau7i;
      } else {
        tau9r = s1 * temp2r - tau7r;
        tau9i = s1 * temp2i - tau7i;
      }
      op[tkm3].re = tau8r - tau9i;
      op[tkm3].im = tau8i + tau9r;
      op[tkm5].re = tau8r + tau9i;
      op[tkm5].im = tau8i - tau9r;
    }
  } else {
    int k, i, ind;
    int M, tkm, u, v, t, tt;
    fft_type temp1r, temp1i, temp2r, temp2i;
    fft_type *wlr = (fft_type *)malloc(sizeof(fft_type) * (radix - 1));
    fft_type *wli = (fft_type *)malloc(sizeof(fft_type) * (radix - 1));
    fft_type *taur = (fft_type *)malloc(sizeof(fft_type) * (radix - 1));
    fft_type *taui = (fft_type *)malloc(sizeof(fft_type) * (radix - 1));
    fft_type *c1 = (fft_type *)malloc(sizeof(fft_type) * (radix - 1));
    fft_type *s1 = (fft_type *)malloc(sizeof(fft_type) * (radix - 1));
    fft_type *yr = (fft_type *)malloc(sizeof(fft_type) * (radix));
    fft_type *yi = (fft_type *)malloc(sizeof(fft_type) * (radix));
    m = N / radix;
    ll = radix * l;
    for (i = 0; i < radix; ++i) {
      mixed_radix_dit_rec(op + i * m, ip + i * l, obj, sgn, m, ll, inc + 1);
    }
    M = (radix - 1) / 2;
    for (i = 1; i < M + 1; ++i) {
      c1[i - 1] = cos(i * PI2 / radix);
      s1[i - 1] = sin(i * PI2 / radix);
    }
    for (i = 0; i < M; ++i) {
      s1[i + M] = -s1[M - 1 - i];
      c1[i + M] = c1[M - 1 - i];
    }
    for (k = 0; k < m; ++k) {
      ind = m - 1 + (radix - 1) * k;
      yr[0] = op[k].re;
      yi[0] = op[k].im;
      for (i = 0; i < radix - 1; ++i) {
        wlr[i] = (obj->twiddle + ind)->re;
        wli[i] = (obj->twiddle + ind)->im;
        tkm = k + (i + 1) * m;
        yr[i + 1] = op[tkm].re * wlr[i] - op[tkm].im * wli[i];
        yi[i + 1] = op[tkm].im * wlr[i] + op[tkm].re * wli[i];
        ind++;
      }
      for (i = 0; i < M; ++i) {
        taur[i] = yr[i + 1] + yr[radix - 1 - i];
        taui[i + M] = yi[i + 1] - yi[radix - 1 - i];
        taui[i] = yi[i + 1] + yi[radix - 1 - i];
        taur[i + M] = yr[i + 1] - yr[radix - 1 - i];
      }
      temp1r = yr[0];
      temp1i = yi[0];
      for (i = 0; i < M; ++i) {
        temp1r += taur[i];
        temp1i += taui[i];
      }
      op[k].re = temp1r;
      op[k].im = temp1i;
      for (u = 0; u < M; u++) {
        temp1r = yr[0];
        temp1i = yi[0];
        temp2r = 0.0;
        temp2i = 0.0;
        for (v = 0; v < M; v++) {
          // int ind2 = (u+v)%M;
          t = (u + 1) * (v + 1);
          while (t >= radix) t -= radix;
          tt = t - 1;
          temp1r += c1[tt] * taur[v];
          temp1i += c1[tt] * taui[v];
          temp2r -= s1[tt] * taur[v + M];
          temp2i -= s1[tt] * taui[v + M];
        }
        temp2r = sgn * temp2r;
        temp2i = sgn * temp2i;
        op[k + (u + 1) * m].re = temp1r - temp2i;
        op[k + (u + 1) * m].im = temp1i + temp2r;
        op[k + (radix - u - 1) * m].re = temp1r + temp2i;
        op[k + (radix - u - 1) * m].im = temp1i - temp2r;
      }
    }
    free(wlr);
    free(wli);
    free(taur);
    free(taui);
    free(c1);
    free(s1);
    free(yr);
    free(yi);
  }
}

void fft_exec(fft_object obj, fft_data *inp, fft_data *oup) {
  if (obj->lt == 0) {
    // fftct_radix3_dit_rec(inp,oup,obj, obj->sgn, obj->N);
    // fftct_mixed_rec(inp,oup,obj, obj->sgn, obj->N);
    // printf("%f \n", 1.785);
    int l, inc;
    int nn, sgn1;
    nn = obj->N;
    sgn1 = obj->sgn;
    l = 1;
    inc = 0;
    // radix3_dit_rec(oup,inp,obj,sgn1,nn,l);
    mixed_radix_dit_rec(oup, inp, obj, sgn1, nn, l, inc);
  } else if (obj->lt == 1) {
    // printf("%f \n", 1.785);
    int nn, sgn1;
    nn = obj->N;
    sgn1 = obj->sgn;
    bluestein_fft(inp, oup, obj, sgn1, nn);
  }
}

void free_fft(fft_object object) { free(object); }

// CWT
cwt_type factorial(int N) {
  // static const cwt_type fact[41] = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320,
  // 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200,
  // 1307674368000,
  //  20922789888000, 355687428096000, 6402373705728000, 121645100408832000,
  //  2432902008176640000, 51090942171709440000.0, 1124000727777607680000.0,
  //  25852016738884976640000.0, 620448401733239439360000.0,
  //  15511210043330985984000000.0, 403291461126605635584000000.0,
  //  10888869450418352160768000000.0, 304888344611713860501504000000.0,
  //  8841761993739701954543616000000.0, 265252859812191058636308480000000.0,
  //  8222838654177922817725562880000000.0,
  //  263130836933693530167218012160000000.0,
  //  8683317618811886495518194401280000000.0,
  //  295232799039604140847618609643520000000.0,
  //  10333147966386144929666651337523200000000.0,
  //  371993326789901217467999448150835200000000.0,
  //  13763753091226345046315979581580902400000000.0,
  //  523022617466601111760007224100074291200000000.0,
  //  20397882081197443358640281739902897356800000000.0,
  //  815915283247897734345611269596115894272000000000.0 };
  static const cwt_type fact[14] = {
      1,    1,     2,      6,       24,       120,       720,
      5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800};
  //  if (N > 40 || N < 0) {
  //    printf("This program is only valid for 0 <= N <= 40 \n");
  //    return -1.0;
  //  }
  if (N > 14 || N < 0) {
    printf("This program is only valid for 0 <= N <= 14 \n");
    return -1.0;
  }
  return fact[N];
}

static cwt_type fix(cwt_type x) {
  // Rounds to the integer nearest to zero
  if (x >= 0.) {
    return floor(x);
  } else {
    return ceil(x);
  }
}

cwt_type cwt_gamma(cwt_type x) {
  /*
   * This C program code is based on  W J Cody's fortran code.
   * http://www.netlib.org/specfun/gamma
   *
   * References:
   "An Overview of Software Development for Special Functions",
   W. J. Cody, Lecture Notes in Mathematics, 506,
   Numerical Analysis Dundee, 1975, G. A. Watson (ed.),
   Springer Verlag, Berlin, 1976.

   Computer Approximations, Hart, Et. Al., Wiley and sons, New York, 1968.
   */

  // numerator and denominator coefficients for 1 <= x <= 2

  cwt_type y, oup, fact, sum, y2, yi, z, nsum, dsum;
  int swi, n, i;
  cwt_type spi = 0.9189385332046727417803297;
  cwt_type pi = 3.1415926535897932384626434;
  cwt_type xmax = 171.624e+0;
  cwt_type xinf = 1.79e308;
  cwt_type eps = 2.22e-16;
  cwt_type xninf = 1.79e-308;
  cwt_type num[8] = {
      -1.71618513886549492533811e+0, 2.47656508055759199108314e+1,
      -3.79804256470945635097577e+2, 6.29331155312818442661052e+2,
      8.66966202790413211295064e+2,  -3.14512729688483675254357e+4,
      -3.61444134186911729807069e+4, 6.64561438202405440627855e+4};
  cwt_type den[8] = {
      -3.08402300119738975254353e+1, 3.15350626979604161529144e+2,
      -1.01515636749021914166146e+3, -3.10777167157231109440444e+3,
      2.25381184209801510330112e+4,  4.75584627752788110767815e+3,
      -1.34659959864969306392456e+5, -1.15132259675553483497211e+5};
  // Coefficients for Hart's Minimax approximation x >= 12
  cwt_type c[7] = {-1.910444077728e-03,
                   8.4171387781295e-04,
                   -5.952379913043012e-04,
                   7.93650793500350248e-04,
                   -2.777777777777681622553e-03,
                   8.333333333333333331554247e-02,
                   5.7083835261e-03};
  y = x;
  swi = 0;
  fact = 1.0;
  n = 0;
  if (y < 0.) {
    // Negative x
    y = -x;
    yi = fix(y);
    oup = y - yi;
    if (oup != 0.0) {
      if (yi != fix(yi * .5) * 2.) {
        swi = 1;
      }
      fact = -pi / sin(pi * oup);
      y += 1.;
    } else {
      return xinf;
    }
  }
  if (y < eps) {
    if (y >= xninf) {
      oup = 1.0 / y;
    } else {
      return xinf;
    }
  } else if (y < 12.) {
    yi = y;
    if (y < 1.) {
      z = y;
      y += 1.;
    } else {
      n = (int)y - 1;
      y -= (cwt_type)n;
      z = y - 1.0;
    }
    nsum = 0.;
    dsum = 1.;
    for (i = 0; i < 8; ++i) {
      nsum = (nsum + num[i]) * z;
      dsum = dsum * z + den[i];
    }
    oup = nsum / dsum + 1.;
    if (yi < y) {
      oup /= yi;
    } else if (yi > y) {
      for (i = 0; i < n; ++i) {
        oup *= y;
        y += 1.;
      }
    }
  } else {
    if (y <= xmax) {
      y2 = y * y;
      sum = c[6];
      for (i = 0; i < 6; ++i) {
        sum = sum / y2 + c[i];
      }
      sum = sum / y - y + spi;
      sum += (y - .5) * log(y);
      oup = exp(sum);
    } else {
      return (xinf);
    }
  }
  if (swi) {
    oup = -oup;
  }
  if (fact != 1.) {
    oup = fact / oup;
  }
  return oup;
}

static void wave_function(int nk, cwt_type dt, int mother, cwt_type param,
                          cwt_type scale1, cwt_type *kwave, cwt_type pi,
                          cwt_type *period1, cwt_type *coi1,
                          fft_data *daughter) {
  cwt_type norm, expnt, fourier_factor;
  int k, m;
  cwt_type temp;
  int sign, re;
  if (mother == 0) {
    // MORLET
    if (param < 0.0) {
      param = 6.0;
    }
    norm = sqrt(2.0 * pi * scale1 / dt) * pow(pi, -0.25);
    for (k = 1; k <= nk / 2 + 1; ++k) {
      temp = (scale1 * kwave[k - 1] - param);
      expnt = -0.5 * temp * temp;
      daughter[k - 1].re = norm * exp(expnt);
      daughter[k - 1].im = 0.0;
    }
    for (k = nk / 2 + 2; k <= nk; ++k) {
      daughter[k - 1].re = daughter[k - 1].im = 0.0;
    }
    fourier_factor = (4.0 * pi) / (param + sqrt(2.0 + param * param));
    *period1 = scale1 * fourier_factor;
    *coi1 = fourier_factor / sqrt(2.0);
  } else if (mother == 1) {
    // PAUL
    if (param < 0.0) {
      param = 4.0;
    }
    m = (int)param;
    norm = sqrt(2.0 * pi * scale1 / dt) *
           (pow(2.0, (cwt_type)m) / sqrt((cwt_type)(m * factorial(2 * m - 1))));
    for (k = 1; k <= nk / 2 + 1; ++k) {
      temp = scale1 * kwave[k - 1];
      expnt = -temp;
      daughter[k - 1].re = norm * pow(temp, (cwt_type)m) * exp(expnt);
      daughter[k - 1].im = 0.0;
    }
    for (k = nk / 2 + 2; k <= nk; ++k) {
      daughter[k - 1].re = daughter[k - 1].im = 0.0;
    }
    fourier_factor = (4.0 * pi) / (2.0 * m + 1.0);
    *period1 = scale1 * fourier_factor;
    *coi1 = fourier_factor * sqrt(2.0);
  } else if (mother == 2) {
    if (param < 0.0) {
      param = 2.0;
    }
    m = (int)param;
    if (m % 2 == 0) {
      re = 1;
    } else {
      re = 0;
    }
    if (m % 4 == 0 || m % 4 == 1) {
      sign = -1;
    } else {
      sign = 1;
    }
    norm = sqrt(2.0 * pi * scale1 / dt) * sqrt(1.0 / cwt_gamma(m + 0.50));
    norm *= sign;
    if (re == 1) {
      for (k = 1; k <= nk; ++k) {
        temp = scale1 * kwave[k - 1];
        daughter[k - 1].re =
            norm * pow(temp, (cwt_type)m) * exp(-0.50 * pow(temp, 2.0));
        daughter[k - 1].im = 0.0;
      }
    } else if (re == 0) {
      for (k = 1; k <= nk; ++k) {
        temp = scale1 * kwave[k - 1];
        daughter[k - 1].re = 0.0;
        daughter[k - 1].im =
            norm * pow(temp, (cwt_type)m) * exp(-0.50 * pow(temp, 2.0));
      }
    }
    fourier_factor = (2.0 * pi) * sqrt(2.0 / (2.0 * m + 1.0));
    *period1 = scale1 * fourier_factor;
    *coi1 = fourier_factor / sqrt(2.0);
  }
}

cwt_object cwt_init(const char *wave, cwt_type param, int siglength,
                    cwt_type dt, int J) {
  cwt_object obj = NULL;
  int N, i, nj2, ibase2, mother = 0;
  cwt_type s0 = 2 * dt, dj;
  cwt_type t1;
  int m, odd;
  const char *pdefault = "pow";
  m = (int)param;
  odd = 1;
  if (2 * (m / 2) == m) {
    odd = 0;
  }
  N = siglength;
  nj2 = 2 * N * J;
  obj = (cwt_object)malloc(sizeof(struct cwt_set) +
                           sizeof(cwt_type) * (nj2 + 2 * J + N));
  if (!strcmp(wave, "morlet") || !strcmp(wave, "morl")) {
    s0 = 2 * dt;
    dj = 0.4875;
    mother = 0;
    if (param < 0.0) {
      printf("\n Morlet Wavelet Parameter should be >= 0 \n");
      exit(-1);
    }
    if (param == 0) {
      param = 6.0;
    }
    strcpy(obj->wave, "morlet");
  } else if (!strcmp(wave, "paul")) {
    s0 = 2 * dt;
    dj = 0.4875;
    mother = 1;
    if (param < 0 || param > 20) {
      printf("\n Paul Wavelet Parameter should be > 0 and <= 20 \n");
      exit(-1);
    }
    if (param == 0) {
      param = 4.0;
    }
    strcpy(obj->wave, "paul");
  } else if (!strcmp(wave, "dgauss") || !strcmp(wave, "dog")) {
    s0 = 2 * dt;
    dj = 0.4875;
    mother = 2;
    if (param < 0 || odd == 1) {
      printf("\n DOG Wavelet Parameter should be > 0 and even \n");
      exit(-1);
    }
    if (param == 0) {
      param = 2.0;
    }
    strcpy(obj->wave, "dog");
  }
  obj->pow = 2;
  strcpy(obj->type, pdefault);
  obj->s0 = s0;
  obj->dj = dj;
  obj->dt = dt;
  obj->J = J;
  obj->siglength = siglength;
  obj->sflag = 0;
  obj->pflag = 1;
  obj->mother = mother;
  obj->m = param;
  t1 = 0.499999 + log((cwt_type)N) / log(2.0);
  ibase2 = 1 + (int)t1;
  obj->npad = (int)pow(2.0, (cwt_type)ibase2);
  obj->output = (cplx_data *)&obj->params[0];
  obj->scale = &obj->params[nj2];
  obj->period = &obj->params[nj2 + J];
  obj->coi = &obj->params[nj2 + 2 * J];
  for (i = 0; i < nj2 + 2 * J + N; ++i) {
    obj->params[i] = 0.0;
  }
  return obj;
}

void setCWTScales(cwt_object wt, cwt_type s0, cwt_type dj, const char *type,
                  int power) {
  int i;
  strcpy(wt->type, type);
  // s0*pow(2.0, (double)(j - 1)*dj);
  if (!strcmp(wt->type, "pow") || !strcmp(wt->type, "power")) {
    for (i = 0; i < wt->J; ++i) {
      wt->scale[i] = s0 * pow((cwt_type)power, (cwt_type)(i)*dj);
    }
    wt->sflag = 1;
    wt->pow = power;
  } else if (!strcmp(wt->type, "lin") || !strcmp(wt->type, "linear")) {
    for (i = 0; i < wt->J; ++i) {
      wt->scale[i] = s0 + (cwt_type)i * dj;
    }
    wt->sflag = 1;
  } else {
    printf("\n Type accepts only two values : pow and lin\n");
    exit(-1);
  }
  wt->s0 = s0;
  wt->dj = dj;
}

void cwt(cwt_object wt, const cwt_type *inp) {
  int i, N, npad, nj2, j, j2;
  N = wt->siglength;
  if (wt->sflag == 0) {
    for (i = 0; i < wt->J; ++i) {
      wt->scale[i] = wt->s0 * pow(2.0, (cwt_type)(i)*wt->dj);
    }
    wt->sflag = 1;
  }
  if (wt->pflag == 0) {
    npad = N;
  } else {
    npad = wt->npad;
  }
  nj2 = 2 * N * wt->J;
  j = wt->J;
  j2 = 2 * j;
  wt->smean = 0.0;
  for (i = 0; i < N; ++i) {
    wt->smean += inp[i];
  }
  wt->smean /= N;
  cwavelet(inp, N, wt->dt, wt->mother, wt->m, wt->s0, wt->dj, wt->J, npad,
           wt->params, wt->params + nj2, wt->params + nj2 + j,
           wt->params + nj2 + j2);
}

void cwavelet(const cwt_type *y, int N, cwt_type dt, int mother, cwt_type param,
              cwt_type s0, cwt_type dj, int jtot, int npad, cwt_type *wave,
              cwt_type *scale, cwt_type *period, cwt_type *coi) {
  int i, j, k, iter;
  cwt_type ymean, freq1, pi, period1, coi1;
  cwt_type tmp1, tmp2;
  cwt_type scale1;
  cwt_type *kwave;
  fft_object obj, iobj;
  fft_data *ypad, *yfft, *daughter;
  (void)s0;
  (void)dj; /* yes, we need these parameters unused */
  pi = 4.0 * atan(1.0);
  if (npad < N) {
    printf("npad must be >= N \n");
    exit(-1);
  }
  obj = fft_init(npad, 1);
  iobj = fft_init(npad, -1);
  ypad = (fft_data *)malloc(sizeof(fft_data) * npad);
  yfft = (fft_data *)malloc(sizeof(fft_data) * npad);
  daughter = (fft_data *)malloc(sizeof(fft_data) * npad);
  kwave = (cwt_type *)malloc(sizeof(cwt_type) * npad);
  ymean = 0.0;
  for (i = 0; i < N; ++i) {
    ymean += y[i];
  }
  ymean /= N;
  for (i = 0; i < N; ++i) {
    ypad[i].re = y[i] - ymean;
    ypad[i].im = 0.0;
  }
  for (i = N; i < npad; ++i) {
    ypad[i].re = ypad[i].im = 0.0;
  }
  // Find FFT of the input y (ypad)
  fft_exec(obj, ypad, yfft);
  for (i = 0; i < npad; ++i) {
    yfft[i].re /= (cwt_type)npad;
    yfft[i].im /= (cwt_type)npad;
  }
  // Construct the wavenumber array
  freq1 = 2.0 * pi / ((cwt_type)npad * dt);
  kwave[0] = 0.0;
  for (i = 1; i < npad / 2 + 1; ++i) {
    kwave[i] = i * freq1;
  }
  for (i = npad / 2 + 1; i < npad; ++i) {
    kwave[i] = -kwave[npad - i];
  }
  // Main loop
  for (j = 1; j <= jtot; ++j) {
    scale1 = scale[j - 1];  // = s0*pow(2.0, (double)(j - 1)*dj);
    wave_function(npad, dt, mother, param, scale1, kwave, pi, &period1, &coi1,
                  daughter);
    period[j - 1] = period1;
    for (k = 0; k < npad; ++k) {
      tmp1 = daughter[k].re * yfft[k].re - daughter[k].im * yfft[k].im;
      tmp2 = daughter[k].re * yfft[k].im + daughter[k].im * yfft[k].re;
      daughter[k].re = tmp1;
      daughter[k].im = tmp2;
    }
    fft_exec(iobj, daughter, ypad);
    iter = 2 * (j - 1) * N;
    for (i = 0; i < N; ++i) {
      wave[iter + 2 * i] = ypad[i].re;
      wave[iter + 2 * i + 1] = ypad[i].im;
    }
  }
  for (i = 1; i <= (N + 1) / 2; ++i) {
    coi[i - 1] = coi1 * dt * ((cwt_type)i - 1.0);
    coi[N - i] = coi[i - 1];
  }
  free(kwave);
  free(ypad);
  free(yfft);
  free(daughter);
  free_fft(obj);
  free_fft(iobj);
}

void cwt_summary(cwt_object wt) {
  printf("\n");
  printf("Wavelet : %s Parameter %lf \n", wt->wave, wt->m);
  printf("\n");
  printf("Length of Input Signal : %d \n", wt->siglength);
  printf("\n");
  printf("Sampling Rate : %g \n", wt->dt);
  printf("\n");
  printf("Total Number of Scales : %d \n", wt->J);
  printf("\n");
  printf("Smallest Scale (s0) : %lf \n", wt->s0);
  printf("\n");
  printf("Separation Between Scales (dj) %lf \n", wt->dj);
  printf("\n");
  printf("Scale Type %s \n", wt->type);
  printf("\n");
  printf(
      "Complex CWT Output Vector is of size %d * %d stored in Row Major format "
      "\n",
      wt->J, wt->siglength);
  printf("\n");
  printf(
      "The ith real value can be accessed using wt->output[i].re and imaginary "
      "value by wt->output[i].im \n");
  printf("\n");
}

void cwt_free(cwt_object object) { free(object); }
}  // namespace sigmap
