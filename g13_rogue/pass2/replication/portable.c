/* NumPy 2.1.2 distribution extracts; see NUMPY_LICENSE.txt and ledger.json.
   Standalone addressed driver v1: no NumPy runtime or global sampler state. */
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#define RAND_INT_TYPE int64_t
#define MIN(x,y) (((x)<(y))?(x):(y))
typedef struct {uint64_t state;} bitgen_t;
static uint64_t mix(uint64_t z) {
  z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);
  z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);
  return z^(z>>31);
}
static double next_double(bitgen_t *s) {
  s->state += UINT64_C(0x9e3779b97f4a7c15);
  return (mix(s->state)>>11)*0x1.0p-53;
}
typedef struct s_binomial_t {
  int has_binomial; /* !=0: following parameters initialized for binomial */
  double psave;
  RAND_INT_TYPE nsave;
  double r;
  double q;
  double fm;
  RAND_INT_TYPE m;
  double p1;
  double xm;
  double xl;
  double xr;
  double c;
  double laml;
  double lamr;
  double p2;
  double p3;
  double p4;
} binomial_t;
double random_loggam(double x) {
  double x0, x2, lg2pi, gl, gl0;
  RAND_INT_TYPE k, n;

  static double a[10] = {8.333333333333333e-02, -2.777777777777778e-03,
                         7.936507936507937e-04, -5.952380952380952e-04,
                         8.417508417508418e-04, -1.917526917526918e-03,
                         6.410256410256410e-03, -2.955065359477124e-02,
                         1.796443723688307e-01, -1.39243221690590e+00};

  if ((x == 1.0) || (x == 2.0)) {
    return 0.0;
  } else if (x < 7.0) {
    n = (RAND_INT_TYPE)(7 - x);
  } else {
    n = 0;
  }
  x0 = x + n;
  x2 = (1.0 / x0) * (1.0 / x0);
  /* log(2 * M_PI) */
  lg2pi = 1.8378770664093453e+00;
  gl0 = a[9];
  for (k = 8; k >= 0; k--) {
    gl0 *= x2;
    gl0 += a[k];
  }
  gl = gl0 / x0 + 0.5 * lg2pi + (x0 - 0.5) * log(x0) - x0;
  if (x < 7.0) {
    for (k = 1; k <= n; k++) {
      gl -= log(x0 - 1.0);
      x0 -= 1.0;
    }
  }
  return gl;
}

static RAND_INT_TYPE random_poisson_mult(bitgen_t *bitgen_state, double lam) {
  RAND_INT_TYPE X;
  double prod, U, enlam;

  enlam = exp(-lam);
  X = 0;
  prod = 1.0;
  while (1) {
    U = next_double(bitgen_state);
    prod *= U;
    if (prod > enlam) {
      X += 1;
    } else {
      return X;
    }
  }
}

static RAND_INT_TYPE random_poisson_ptrs(bitgen_t *bitgen_state, double lam) {
  RAND_INT_TYPE k;
  double U, V, slam, loglam, a, b, invalpha, vr, us;

  slam = sqrt(lam);
  loglam = log(lam);
  b = 0.931 + 2.53 * slam;
  a = -0.059 + 0.02483 * b;
  invalpha = 1.1239 + 1.1328 / (b - 3.4);
  vr = 0.9277 - 3.6224 / (b - 2);

  while (1) {
    U = next_double(bitgen_state) - 0.5;
    V = next_double(bitgen_state);
    us = 0.5 - fabs(U);
    k = (RAND_INT_TYPE)floor((2 * a / us + b) * U + lam + 0.43);
    if ((us >= 0.07) && (V <= vr)) {
      return k;
    }
    if ((k < 0) || ((us < 0.013) && (V > us))) {
      continue;
    }
    /* log(V) == log(0.0) ok here */
    /* if U==0.0 so that us==0.0, log is ok since always returns */
    if ((log(V) + log(invalpha) - log(a / (us * us) + b)) <=
        (-lam + k * loglam - random_loggam(k + 1))) {
      return k;
    }
  }
}

RAND_INT_TYPE random_poisson(bitgen_t *bitgen_state, double lam) {
  if (lam >= 10) {
    return random_poisson_ptrs(bitgen_state, lam);
  } else if (lam == 0) {
    return 0;
  } else {
    return random_poisson_mult(bitgen_state, lam);
  }
}

RAND_INT_TYPE random_binomial_btpe(bitgen_t *bitgen_state, RAND_INT_TYPE n,
                                   double p, binomial_t *binomial) {
  double r, q, fm, p1, xm, xl, xr, c, laml, lamr, p2, p3, p4;
  double a, u, v, s, F, rho, t, A, nrq, x1, x2, f1, f2, z, z2, w, w2, x;
  RAND_INT_TYPE m, y, k, i;

  if (!(binomial->has_binomial) || (binomial->nsave != n) ||
      (binomial->psave != p)) {
    /* initialize */
    binomial->nsave = n;
    binomial->psave = p;
    binomial->has_binomial = 1;
    binomial->r = r = MIN(p, 1.0 - p);
    binomial->q = q = 1.0 - r;
    binomial->fm = fm = n * r + r;
    binomial->m = m = (RAND_INT_TYPE)floor(binomial->fm);
    binomial->p1 = p1 = floor(2.195 * sqrt(n * r * q) - 4.6 * q) + 0.5;
    binomial->xm = xm = m + 0.5;
    binomial->xl = xl = xm - p1;
    binomial->xr = xr = xm + p1;
    binomial->c = c = 0.134 + 20.5 / (15.3 + m);
    a = (fm - xl) / (fm - xl * r);
    binomial->laml = laml = a * (1.0 + a / 2.0);
    a = (xr - fm) / (xr * q);
    binomial->lamr = lamr = a * (1.0 + a / 2.0);
    binomial->p2 = p2 = p1 * (1.0 + 2.0 * c);
    binomial->p3 = p3 = p2 + c / laml;
    binomial->p4 = p4 = p3 + c / lamr;
  } else {
    r = binomial->r;
    q = binomial->q;
    fm = binomial->fm;
    m = binomial->m;
    p1 = binomial->p1;
    xm = binomial->xm;
    xl = binomial->xl;
    xr = binomial->xr;
    c = binomial->c;
    laml = binomial->laml;
    lamr = binomial->lamr;
    p2 = binomial->p2;
    p3 = binomial->p3;
    p4 = binomial->p4;
  }

/* sigh ... */
Step10:
  nrq = n * r * q;
  u = next_double(bitgen_state) * p4;
  v = next_double(bitgen_state);
  if (u > p1)
    goto Step20;
  y = (RAND_INT_TYPE)floor(xm - p1 * v + u);
  goto Step60;

Step20:
  if (u > p2)
    goto Step30;
  x = xl + (u - p1) / c;
  v = v * c + 1.0 - fabs(m - x + 0.5) / p1;
  if (v > 1.0)
    goto Step10;
  y = (RAND_INT_TYPE)floor(x);
  goto Step50;

Step30:
  if (u > p3)
    goto Step40;
  y = (RAND_INT_TYPE)floor(xl + log(v) / laml);
  /* Reject if v==0.0 since previous cast is undefined */
  if ((y < 0) || (v == 0.0))
    goto Step10;
  v = v * (u - p2) * laml;
  goto Step50;

Step40:
  y = (RAND_INT_TYPE)floor(xr - log(v) / lamr);
  /* Reject if v==0.0 since previous cast is undefined */
  if ((y > n) || (v == 0.0))
    goto Step10;
  v = v * (u - p3) * lamr;

Step50:
  k = llabs(y - m);
  if ((k > 20) && (k < ((nrq) / 2.0 - 1)))
    goto Step52;

  s = r / q;
  a = s * (n + 1);
  F = 1.0;
  if (m < y) {
    for (i = m + 1; i <= y; i++) {
      F *= (a / i - s);
    }
  } else if (m > y) {
    for (i = y + 1; i <= m; i++) {
      F /= (a / i - s);
    }
  }
  if (v > F)
    goto Step10;
  goto Step60;

Step52:
  rho =
      (k / (nrq)) * ((k * (k / 3.0 + 0.625) + 0.16666666666666666) / nrq + 0.5);
  t = -k * k / (2 * nrq);
  /* log(0.0) ok here */
  A = log(v);
  if (A < (t - rho))
    goto Step60;
  if (A > (t + rho))
    goto Step10;

  x1 = y + 1;
  f1 = m + 1;
  z = n + 1 - m;
  w = n - y + 1;
  x2 = x1 * x1;
  f2 = f1 * f1;
  z2 = z * z;
  w2 = w * w;
  if (A > (xm * log(f1 / x1) + (n - m + 0.5) * log(z / w) +
           (y - m) * log(w * r / (x1 * q)) +
           (13680. - (462. - (132. - (99. - 140. / f2) / f2) / f2) / f2) / f1 /
               166320. +
           (13680. - (462. - (132. - (99. - 140. / z2) / z2) / z2) / z2) / z /
               166320. +
           (13680. - (462. - (132. - (99. - 140. / x2) / x2) / x2) / x2) / x1 /
               166320. +
           (13680. - (462. - (132. - (99. - 140. / w2) / w2) / w2) / w2) / w /
               166320.)) {
    goto Step10;
  }

Step60:
  if (p > 0.5) {
    y = n - y;
  }

  return y;
}

RAND_INT_TYPE random_binomial_inversion(bitgen_t *bitgen_state, RAND_INT_TYPE n,
                                        double p, binomial_t *binomial) {
  double q, qn, np, px, U;
  RAND_INT_TYPE X, bound;

  if (!(binomial->has_binomial) || (binomial->nsave != n) ||
      (binomial->psave != p)) {
    binomial->nsave = n;
    binomial->psave = p;
    binomial->has_binomial = 1;
    binomial->q = q = 1.0 - p;
    binomial->r = qn = exp(n * log(q));
    binomial->c = np = n * p;
    binomial->m = bound = (RAND_INT_TYPE)MIN(n, np + 10.0 * sqrt(np * q + 1));
  } else {
    q = binomial->q;
    qn = binomial->r;
    np = binomial->c;
    bound = binomial->m;
  }
  X = 0;
  px = qn;
  U = next_double(bitgen_state);
  while (U > px) {
    X++;
    if (X > bound) {
      X = 0;
      px = qn;
      U = next_double(bitgen_state);
    } else {
      U -= px;
      px = ((n - X + 1) * p * px) / (X * q);
    }
  }
  return X;
}

int64_t random_binomial(bitgen_t *bitgen_state, double p, int64_t n,
                        binomial_t *binomial) {
  double q;

  if ((n == 0LL) || (p == 0.0f))
    return 0;

  if (p <= 0.5) {
    if (p * n <= 30.0) {
      return random_binomial_inversion(bitgen_state, n, p, binomial);
    } else {
      return random_binomial_btpe(bitgen_state, n, p, binomial);
    }
  } else {
    q = 1.0 - p;
    if (q * n <= 30.0) {
      return n - random_binomial_inversion(bitgen_state, n, q, binomial);
    } else {
      return n - random_binomial_btpe(bitgen_state, n, q, binomial);
    }
  }
}

/* Each scalar distribution has a fresh addressed rejection substream. */
void portable_binomial(int64_t size, const int64_t *n, const double *p,
                       const uint64_t *keys, int64_t *out) {
  for(int64_t i=0;i<size;i++) {
    bitgen_t s={keys[i]}; binomial_t b={0};
    double q=floor(p[i]*0x1.0p40+0.5)*0x1.0p-40;
    out[i]=random_binomial(&s,q,n[i],&b);
  }
}
int64_t portable_poisson(double lam, uint64_t key) {
  bitgen_t s={key}; return random_poisson(&s,lam);
}
void portable_uniform(int64_t size, const uint64_t *keys, double *out) {
  for(int64_t i=0;i<size;i++) {bitgen_t s={keys[i]}; out[i]=next_double(&s);}
}
void portable_multinomial(int64_t n, int64_t size, const double *p,
                           uint64_t key, int64_t *out) {
  double remaining=1.;
  for(int64_t i=0;i<size-1;i++) {
    bitgen_t s={mix(key ^ mix((uint64_t)i))}; binomial_t b={0};
    /* The sequential conditional probabilities are quantized, not the counts. */
    double q=floor(fmin(1.,fmax(0.,p[i]/remaining))*0x1.0p40+0.5)*0x1.0p-40;
    out[i]=random_binomial(&s,q,n,&b); n-=out[i]; remaining-=p[i];
  }
  out[size-1]=n;
}
