//#include <fem.hpp> // Fortran EMulation library of fable module
#include <string>

namespace opendap_convert {

using namespace major_types;

void
ahouse2(...)
{
  throw std::runtime_error(
    "Missing function implementation: ahouse2");
}

void
normylm(...)
{
  throw std::runtime_error(
    "Missing function implementation: normylm");
}

void
ylm(...)
{
  throw std::runtime_error(
    "Missing function implementation: ylm");
}

struct common_splhprm
{
  static const int mxknt = 21;

  float spknt[];
  arr<float, 2> qq0;
  arr<float, 3> qq;

  common_splhprm() :
    spknt(dimension(mxknt), fem::fill0),
    qq0(dimension(mxknt, mxknt), fem::fill0),
    qq(dimension(3, mxknt, mxknt), fem::fill0)
  {}
};

const int common_splhprm::mxknt;

struct common :
  common,
  common_splhprm
{
  fem::cmn_sve rspln_sve;
  fem::cmn_sve splhsetup_sve;
  fem::cmn_sve rsple_sve;

  common(
    int argc,
    char const* argv[])
  :
    fem::common(argc, argv)
  {}
};

//C
//C ------------------------------------------------------------------------------
//C
int
istlen(
  str_cref string)
{
  int return_value = fem::int0;
  fem::str<1> null = fem::char0;
  int k = fem::int0;
  int i = fem::int0;
  int j = fem::int0;
  null = fem::fchar(0);
  k = fem::len(string);
  FEM_DO_SAFE(i, 1, k) {
    j = k + 1 - i;
    if (string(j, j) == " " || string(j, j) == null) {
      goto statement_10;
    }
    return_value = j;
    goto statement_99;
    statement_10:;
  }
  return_value = 0;
  statement_99:
  return return_value;
}

//C
//C ---------------------------------------------------------------------
//C
void
wsphhead(
  common& cmn,
  int const& lmax,
  int const& idp1,
  int const& idp2,
  str_ref strt,
  int& lstrt)
{
  common_write write(cmn);
  //C
  const int maxl = 100;
  if (lmax > maxl) {
    FEM_STOP("wsphhead >>> lmax too big for stment");
  }
  const int maxp = 24;
  if (idp2 > maxp) {
    FEM_STOP("wsphhead >>> idp2 too large");
  }
  int i = fem::int0;
  arr_1d<101, int> lask(dim1(0, maxl), fem::fill0);
  FEM_DO_SAFE(i, 0, maxl) {
    if (i <= lmax) {
      lask(i) = 1;
    }
    else {
      lask(i) = 0;
    }
  }
  //C
  arr_1d<maxp, int> mask(fem::fill0);
  FEM_DO_SAFE(i, 1, 24) {
    mask(i) = 0;
  }
  FEM_DO_SAFE(i, idp1, idp2) {
    mask(i) = 1;
  }
  //C
  fem::str<120> str1 = fem::char0;
  {
    write_loop wloop(str1, "(10x,i4,1x,100i1)");
    wloop, lmax;
    FEM_DO_SAFE(i, 0, lmax) {
      wloop, lask(i);
    }
  }
  fem::str<120> str2 = fem::char0;
  {
    write_loop wloop(str2, "(i4,1x,24i1)");
    wloop, 24;
    FEM_DO_SAFE(i, 1, 24) {
      wloop, mask(i);
    }
  }
  //C
  int len1 = istlen(str1);
  int len2 = istlen(str2);
  if ((len1 + len2) > 120) {
    FEM_STOP("wsphhead >>> increase dim of strt");
  }
  strt = str1(1, len1) + str2(1, len2);
  lstrt = len1 + len2 + 1;
}

struct rspln_save
{
  fem::variant_core_and_bindings save_equivalences;
  float small;

  rspln_save() :
    small(fem::float0)
  {}
};

//C
//C --------------------------------------------------------
//C
void
rspln(
  common& cmn,
  int const& i1,
  int const& i2,
  arr_cref<float> x,
  arr_cref<float> y,
  arr_ref<float, 2> q,
  arr_ref<float, 2> f)
{
  FEM_CMN_SVE(rspln);
  x(dimension(star));
  y(dimension(star));
  q(dimension(3, star));
  f(dimension(3, star));
  save_equivalences sve_equivalences(sve.save_equivalences);
  float& small = sve.small;
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<float> yy(dimension(3));
      mbr<float> y0;
      sve_equivalences.allocate(),
        equivalence(yy, y0)
          .align<1>(arr_index(1))
           .with<2>()
      ;
    }
  }
  arr_ref<float> yy(sve_equivalences.bind<float>(), dimension(3));
  float& y0 = sve_equivalences.bind<float>();
  if (is_called_first_time) {
    small = 1.e-5f;
    {
      static const float values[] = {
        0.f, 0.f, 0.f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        yy;
    }
  }
  int j1 = fem::int0;
  float a0 = fem::float0;
  int i = fem::int0;
  float b0 = fem::float0;
  int j2 = fem::int0;
  int j = fem::int0;
  float h = fem::float0;
  float h2 = fem::float0;
  float b1 = fem::float0;
  float ha = fem::float0;
  float h2a = fem::float0;
  float h3a = fem::float0;
  float h2b = fem::float0;
  int k = fem::int0;
  //C
  //C C$C$C$C$C$ CALLS ONLY LIBRARY ROUTINES C$C$C$C$C$
  //C
  //C   SUBROUTINE RSPLN COMPUTES CUBIC SPLINE INTERPOLATION COEFFICIENTS
  //C   FOR Y(X) BETWEEN GRID POINTS I1 AND I2 SAVING THEM IN Q.  THE
  //C   INTERPOLATION IS CONTINUOUS WITH CONTINUOUS FIRST AND SECOND
  //C   DERIVITIVES.  IT AGREES EXACTLY WITH Y AT GRID POINTS AND WITH THE
  //C   THREE POINT FIRST DERIVITIVES AT BOTH END POINTS (I1 AND I2).
  //C   X MUST BE MONOTONIC BUT IF TWO SUCCESSIVE VALUES OF X ARE EQUAL
  //C   A DISCONTINUITY IS ASSUMED AND SEPERATE INTERPOLATION IS DONE ON
  //C   EACH STRICTLY MONOTONIC SEGMENT.  THE ARRAYS MUST BE DIMENSIONED AT
  //C   LEAST - X(I2), Y(I2), Q(3,I2), AND F(3,I2).  F IS WORKING STORAGE
  //C   FOR RSPLN.
  //C                                                     -RPB
  //C
  j1 = i1 + 1;
  y0 = 0.f;
  //C   BAIL OUT IF THERE ARE LESS THAN TWO POINTS TOTAL.
  switch (fem::if_arithmetic(i2 - i1)) {
    case -1: goto statement_13;
    case  0: goto statement_17;
    default: goto statement_8;
  }
  statement_8:
  a0 = x(j1 - 1);
  //C   SEARCH FOR DISCONTINUITIES.
  FEM_DO_SAFE(i, j1, i2) {
    b0 = a0;
    a0 = x(i);
    if (fem::abs((a0 - b0) / fem::amax1(a0, b0)) < small) {
      goto statement_4;
    }
  }
  statement_17:
  j1 = j1 - 1;
  j2 = i2 - 2;
  goto statement_5;
  statement_4:
  j1 = j1 - 1;
  j2 = i - 3;
  //C   SEE IF THERE ARE ENOUGH POINTS TO INTERPOLATE (AT LEAST THREE).
  statement_5:
  switch (fem::if_arithmetic(j2 + 1 - j1)) {
    case -1: goto statement_9;
    case  0: goto statement_10;
    default: goto statement_11;
  }
  //C   ONLY TWO POINTS.  USE LINEAR INTERPOLATION.
  statement_10:
  j2 += 2;
  y0 = (y(j2) - y(j1)) / (x(j2) - x(j1));
  FEM_DO_SAFE(j, 1, 3) {
    q(j, j1) = yy(j);
    q(j, j2) = yy(j);
  }
  goto statement_12;
  //C   MORE THAN TWO POINTS.  DO SPLINE INTERPOLATION.
  statement_11:
  a0 = 0.f;
  h = x(j1 + 1) - x(j1);
  h2 = x(j1 + 2) - x(j1);
  y0 = h * h2 * (h2 - h);
  h = h * h;
  h2 = h2 * h2;
  //C   CALCULATE DERIVITIVE AT NEAR END.
  b0 = (y(j1) * (h - h2) + y(j1 + 1) * h2 - y(j1 + 2) * h) / y0;
  b1 = b0;
  //C   EXPLICITLY REDUCE BANDED MATRIX TO AN UPPER BANDED MATRIX.
  FEM_DO_SAFE(i, j1, j2) {
    h = x(i + 1) - x(i);
    y0 = y(i + 1) - y(i);
    h2 = h * h;
    ha = h - a0;
    h2a = h - 2.f * a0;
    h3a = 2.f * h - 3.f * a0;
    h2b = h2 * b0;
    q(1, i) = h2 / ha;
    q(2, i) = -ha / (h2a * h2);
    q(3, i) = -h * h2a / h3a;
    f(1, i) = (y0 - h * b0) / (h * ha);
    f(2, i) = (h2b - y0 * (2.f * h - a0)) / (h * h2 * h2a);
    f(3, i) = -(h2b - 3.f * y0 * ha) / (h * h3a);
    a0 = q(3, i);
    b0 = f(3, i);
  }
  //C   TAKE CARE OF LAST TWO ROWS.
  i = j2 + 1;
  h = x(i + 1) - x(i);
  y0 = y(i + 1) - y(i);
  h2 = h * h;
  ha = h - a0;
  h2a = h * ha;
  h2b = h2 * b0 - y0 * (2.f * h - a0);
  q(1, i) = h2 / ha;
  f(1, i) = (y0 - h * b0) / h2a;
  ha = x(j2) - x(i + 1);
  y0 = -h * ha * (ha + h);
  ha = ha * ha;
  //C   CALCULATE DERIVITIVE AT FAR END.
  y0 = (y(i + 1) * (h2 - ha) + y(i) * ha - y(j2) * h2) / y0;
  q(3, i) = (y0 * h2a + h2b) / (h * h2 * (h - 2.f * a0));
  q(2, i) = f(1, i) - q(1, i) * q(3, i);
  //C   SOLVE UPPER BANDED MATRIX BY REVERSE ITERATION.
  FEM_DO_SAFE(j, j1, j2) {
    k = i - 1;
    q(1, i) = f(3, k) - q(3, k) * q(2, i);
    q(3, k) = f(2, k) - q(2, k) * q(1, i);
    q(2, k) = f(1, k) - q(1, k) * q(3, k);
    i = k;
  }
  q(1, i) = b1;
  //C   FILL IN THE LAST POINT WITH A LINEAR EXTRAPOLATION.
  statement_9:
  j2 += 2;
  FEM_DO_SAFE(j, 1, 3) {
    q(j, j2) = yy(j);
  }
  //C   SEE IF THIS DISCONTINUITY IS THE LAST.
  statement_12:
  switch (fem::if_arithmetic(j2 - i2)) {
    case -1: goto statement_6;
    case  0: goto statement_13;
    default: goto statement_13;
  }
  //C   NO.  GO BACK FOR MORE.
  statement_6:
  j1 = j2 + 2;
  switch (fem::if_arithmetic(j1 - i2)) {
    case -1: goto statement_8;
    case  0: goto statement_8;
    default: goto statement_7;
  }
  //C   THERE IS ONLY ONE POINT LEFT AFTER THE LATEST DISCONTINUITY.
  statement_7:
  FEM_DO_SAFE(j, 1, 3) {
    q(j, i2) = yy(j);
  }
  //C   FINI.
  statement_13:;
}

struct splhsetup_save
{
};

//C
//C ------------------------------------------------------------------------------
//C
void
splhsetup(
  common& cmn)
{
  FEM_CMN_SVE(splhsetup);
  // COMMON splhprm
  const int mxknt = 21;
  arr_ref<float> spknt(cmn.spknt, dimension(mxknt));
  arr_ref<float, 2> qq0(cmn.qq0, dimension(mxknt, mxknt));
  arr_ref<float, 3> qq(cmn.qq, dimension(3, mxknt, mxknt));
  //
  if (is_called_first_time) {
    static const float values[] = {
      -1.00000f, -0.78631f, -0.59207f, -0.41550f, -0.25499f,
        -0.10909f, 0.02353f, 0.14409f, 0.25367f, 0.35329f, 0.44384f,
        0.52615f, 0.60097f, 0.66899f, 0.73081f, 0.78701f, 0.83810f,
        0.88454f, 0.92675f, 0.96512f, 1.00000f
    };
    fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
      spknt;
  }
  int i = fem::int0;
  int j = fem::int0;
  FEM_DO_SAFE(i, 1, mxknt) {
    FEM_DO_SAFE(j, 1, mxknt) {
      if (i == j) {
        qq0(j, i) = 1.f;
      }
      else {
        qq0(j, i) = 0.f;
      }
    }
  }
  arr_2d<3, mxknt, float> qqwk(fem::fill0);
  FEM_DO_SAFE(i, 1, mxknt) {
    rspln(cmn, 1, mxknt, spknt(1), qq0(1, i), qq(1, 1, i), qqwk(1, 1));
  }
}

struct rsple_save
{
  int i;

  rsple_save() :
    i(fem::int0)
  {}
};

//C
//C -------------------------------------------------
//C
float
rsple(
  common& cmn,
  int const& i1,
  int const& i2,
  arr_cref<float> x,
  arr_cref<float> y,
  arr_cref<float, 2> q,
  float const& s)
{
  float return_value = fem::float0;
  FEM_CMN_SVE(rsple);
  x(dimension(star));
  y(dimension(star));
  q(dimension(3, star));
  int& i = sve.i;
  if (is_called_first_time) {
    i = 1;
  }
  int ii = fem::int0;
  float h = fem::float0;
  //C
  //C C$C$C$C$C$ CALLS ONLY LIBRARY ROUTINES C$C$C$C$C$
  //C
  //C   RSPLE RETURNS THE VALUE OF THE FUNCTION Y(X) EVALUATED AT POINT S
  //C   USING THE CUBIC SPLINE COEFFICIENTS COMPUTED BY RSPLN AND SAVED IN
  //C   Q.  IF S IS OUTSIDE THE INTERVAL (X(I1),X(I2)) RSPLE EXTRAPOLATES
  //C   USING THE FIRST OR LAST INTERPOLATION POLYNOMIAL.  THE ARRAYS MUST
  //C   BE DIMENSIONED AT LEAST - X(I2), Y(I2), AND Q(3,I2).
  //C
  //C                                                     -RPB
  ii = i2 - 1;
  //C   GUARANTEE I WITHIN BOUNDS.
  i = fem::max0(i, i1);
  i = fem::min0(i, ii);
  //C   SEE IF X IS INCREASING OR DECREASING.
  switch (fem::if_arithmetic(x(i2) - x(i1))) {
    case -1: goto statement_1;
    case  0: goto statement_2;
    default: goto statement_2;
  }
  //C   X IS DECREASING.  CHANGE I AS NECESSARY.
  statement_1:
  switch (fem::if_arithmetic(s - x(i))) {
    case -1: goto statement_3;
    case  0: goto statement_3;
    default: goto statement_4;
  }
  statement_4:
  i = i - 1;
  switch (fem::if_arithmetic(i - i1)) {
    case -1: goto statement_11;
    case  0: goto statement_6;
    default: goto statement_1;
  }
  statement_3:
  switch (fem::if_arithmetic(s - x(i + 1))) {
    case -1: goto statement_5;
    case  0: goto statement_6;
    default: goto statement_6;
  }
  statement_5:
  i++;
  switch (fem::if_arithmetic(i - ii)) {
    case -1: goto statement_3;
    case  0: goto statement_6;
    default: goto statement_7;
  }
  //C   X IS INCREASING.  CHANGE I AS NECESSARY.
  statement_2:
  switch (fem::if_arithmetic(s - x(i + 1))) {
    case -1: goto statement_8;
    case  0: goto statement_8;
    default: goto statement_9;
  }
  statement_9:
  i++;
  switch (fem::if_arithmetic(i - ii)) {
    case -1: goto statement_2;
    case  0: goto statement_6;
    default: goto statement_7;
  }
  statement_8:
  switch (fem::if_arithmetic(s - x(i))) {
    case -1: goto statement_10;
    case  0: goto statement_6;
    default: goto statement_6;
  }
  statement_10:
  i = i - 1;
  switch (fem::if_arithmetic(i - i1)) {
    case -1: goto statement_11;
    case  0: goto statement_6;
    default: goto statement_8;
  }
  statement_7:
  i = ii;
  goto statement_6;
  statement_11:
  i = i1;
  //C   CALCULATE RSPLE USING SPLINE COEFFICIENTS IN Y AND Q.
  statement_6:
  h = s - x(i);
  return_value = y(i) + h * (q(1, i) + h * (q(2, i) + h * q(3, i)));
  return return_value;
}

//C
//C ------------------------------------------------
float
splh(
  common& cmn,
  int const& ind,
  float const& x)
{
  float return_value = fem::float0;
  const int mxknt = 21;
  arr_cref<float> spknt(cmn.spknt, dimension(mxknt));
  arr_cref<float, 2> qq0(cmn.qq0, dimension(mxknt, mxknt));
  arr_cref<float, 3> qq(cmn.qq, dimension(3, mxknt, mxknt));
  //
  //C
  if (x > 1 || x <  - 1) {
    return_value = 0.f;
    goto statement_10;
  }
  //C
  return_value = rsple(cmn, 1, mxknt, spknt(1), qq0(1, mxknt - ind),
    qq(1, 1, mxknt - ind), x);
  //C
  statement_10:
  return return_value;
}

//C
//C ------------------------------------------------------------------------------
//C
double
dot(
  arr_cref<double, 2> a,
  int const& j,
  arr_cref<double, 2> b,
  int const& k,
  int const& n)
{
  double return_value = fem::double0;
  a(dimension(j, n));
  b(dimension(k, n));
  return_value = 0.e0;
  int i = fem::int0;
  FEM_DO_SAFE(i, 1, n) {
    return_value += a(1, i) * b(1, i);
  }
  return return_value;
}

//C
//C ------------------------------------------------------------
//C
void
ahouse(
  common& cmn,
  int const& n,
  int const& iu,
  arr_ref<double> c,
  arr_ref<double> y,
  arr_ref<double> a,
  arr_ref<double> b,
  arr_ref<double> p,
  arr_ref<double> ta,
  arr_ref<double> tb,
  arr_ref<double> w,
  arr_ref<double> v,
  arr_ref<double> ev)
{
  c(dimension(star));
  y(dimension(star));
  a(dimension(star));
  b(dimension(star));
  p(dimension(star));
  ta(dimension(star));
  tb(dimension(star));
  w(dimension(star));
  v(dimension(star));
  ev(dimension(star));
  common_write write(cmn);
  double eps = fem::double0;
  double umeps = fem::double0;
  double tol = fem::double0;
  int jskip = fem::int0;
  int kskip = fem::int0;
  int nm1 = fem::int0;
  int i = fem::int0;
  int idm1 = fem::int0;
  double ckj = fem::double0;
  int ip1 = fem::int0;
  int nmi = fem::int0;
  int kj = fem::int0;
  int j = fem::int0;
  int jp1 = fem::int0;
  double vj = fem::double0;
  int k = fem::int0;
  int lj = fem::int0;
  int jd = fem::int0;
  double pj = fem::double0;
  double wj = fem::double0;
  double dc = fem::double0;
  double sum = fem::double0;
  double s = fem::double0;
  double csd = fem::double0;
  double h = fem::double0;
  double sp = fem::double0;
  double u = fem::double0;
  double bd = fem::double0;
  double rbd = fem::double0;
  int ik = fem::int0;
  int ndim = fem::int0;
  double el = fem::double0;
  double elam = fem::double0;
  double du = fem::double0;
  int iag = fem::int0;
  double q = fem::double0;
  double x = fem::double0;
  int m = fem::int0;
  int mm = fem::int0;
  int ii = fem::int0;
  int l = fem::int0;
  double f = fem::double0;
  double t = fem::double0;
  double ay = fem::double0;
  int id = fem::int0;
  double xnorm = fem::double0;
  //C   SCM VERSION
  eps = 1.e-14;
  umeps = 1.e0 - eps;
  tol = 1.e-70;
  jskip = 0;
  kskip = 1;
  nm1 = n - 1;
  i = 1;
  idm1 = 0;
  p(1) = 0.e0;
  v(1) = 0.e0;
  w(1) = 0.e0;
  if (n <= 0) {
    return;
  }
  if (n > 2) {
    goto statement_4;
  }
  if (n == 2) {
    goto statement_3;
  }
  ev(1) = c(1);
  y(1) = 1.e0;
  write(iu, fem::unformatted), ev(1), y(1);
  cmn.io.endfile(iu);
  cmn.io.rewind(iu);
  return;
  statement_3:
  a(1) = c(1);
  b(1) = -c(2);
  ckj = c(3);
  ip1 = 2;
  goto statement_215;
  statement_4:
  ip1 = i + 1;
  nmi = n - i;
  kj = idm1;
  j = i;
  statement_5:
  jp1 = j + 1;
  vj = v(j);
  k = j;
  lj = n - j + 1;
  jd = kj + 1;
  if (kskip == 1) {
    goto statement_6;
  }
  pj = p(j);
  wj = w(j);
  statement_6:
  kj++;
  ckj = c(kj);
  if (kskip == 1) {
    goto statement_7;
  }
  dc = -(pj * w(k) + wj * p(k));
  ckj += dc;
  c(kj) = ckj;
  statement_7:
  if (j > i) {
    goto statement_14;
  }
  if (k > j) {
    goto statement_8;
  }
  a(i) = ckj;
  k++;
  goto statement_6;
  statement_8:
  y(k) = 0.e0;
  v(k) = ckj;
  k++;
  if (k <= n) {
    goto statement_6;
  }
  jskip = 0;
  sum = dot(v(jp1), 1, v(jp1), 1, lj - 1);
  if (sum <= tol) {
    goto statement_10;
  }
  s = fem::dsqrt(sum);
  csd = v(jp1);
  if (csd < 0.e0) {
    s = -s;
  }
  v(jp1) = csd + s;
  c(jd + 1) = v(jp1);
  h = sum + csd * s;
  b(i) = -s;
  goto statement_12;
  statement_10:
  b(i) = 0.e0;
  jskip = 1;
  statement_12:
  idm1 = kj;
  if (jskip == 1 && kskip == 1) {
    goto statement_215;
  }
  j = jp1;
  goto statement_5;
  statement_14:
  if (jskip == 0) {
    goto statement_15;
  }
  k++;
  if (k <= n) {
    goto statement_6;
  }
  j = jp1;
  if (j <= n) {
    goto statement_5;
  }
  goto statement_215;
  statement_15:
  y(k) += ckj * vj;
  k++;
  if (k <= n) {
    goto statement_6;
  }
  if (j == n) {
    goto statement_17;
  }
  y(j) += dot(c(jd + 1), 1, v(jp1), 1, lj - 1);
  j = jp1;
  goto statement_5;
  //C
  statement_17:
  sp = dot(v(ip1), 1, y(ip1), 1, nmi) / (h + h);
  FEM_DO_SAFE(j, ip1, n) {
    w(j) = v(j);
    p(j) = (y(j) - sp * v(j)) / h;
  }
  statement_215:
  kskip = jskip;
  i = ip1;
  if (i <= nm1) {
    goto statement_4;
  }
  a(n) = ckj;
  b(nm1) = -b(nm1);
  b(n) = 0.e0;
  u = fem::dabs(a(1)) + fem::dabs(b(1));
  FEM_DO_SAFE(i, 2, n) {
    u = fem::dmax1(u, fem::dabs(a(i)) + fem::dabs(b(i)) + fem::dabs(b(i - 1)));
  }
  bd = u;
  rbd = 1.e0 / u;
  FEM_DO_SAFE(i, 1, n) {
    w(i) = b(i);
    b(i) = fem::pow2((b(i) / u));
    a(i) = a(i) / u;
    v(i) = 0.e0;
    ev(i) = -1.e0;
  }
  u = 1.e0;
  ik = 1;
  ndim = kj;
  statement_1000:
  k = ik;
  el = ev(k);
  statement_24:
  elam = .5e0 * (u + el);
  du = (4.e0 * fem::dabs(elam) + rbd) * eps;
  if (fem::dabs(u - el) <= du) {
    goto statement_42;
  }
  iag = 0;
  q = a(1) - elam;
  if (q >= 0.e0) {
    iag++;
  }
  FEM_DO_SAFE(i, 2, n) {
    if (q == 0.e0) {
      x = fem::dabs(w(i - 1) / bd) / eps;
    }
    if (q != 0.e0) {
      x = b(i - 1) / q;
    }
    q = a(i) - elam - x;
    if (q >= 0.e0) {
      iag++;
    }
  }
  if (iag >= k) {
    goto statement_39;
  }
  u = elam;
  goto statement_24;
  statement_39:
  if (iag == k) {
    goto statement_41;
  }
  m = k + 1;
  FEM_DO_SAFE(mm, m, iag) {
    ev(mm) = elam;
  }
  statement_41:
  el = elam;
  goto statement_24;
  statement_42:
  elam = bd * elam;
  ev(k) = elam;
  if (ik == 1) {
    goto statement_44;
  }
  if (elam >= ev(ik - 1)) {
    ev(ik) = umeps * ev(ik - 1);
  }
  statement_44:
  i = ik;
  ii = 1;
  l = n - 1;
  FEM_DO_SAFE(j, 1, n) {
    y(j) = 1.e0;
  }
  statement_50:
  FEM_DO_SAFE(k, 1, n) {
    p(k) = 0.e0;
    tb(k) = w(k);
    ta(k) = bd * a(k) - ev(i);
  }
  j = 1;
  FEM_DO_SAFE(jp1, 2, n) {
    if (fem::dabs(ta(j)) < fem::dabs(w(j))) {
      goto statement_53;
    }
    if (ta(j) == 0.e0) {
      ta(j) = eps;
    }
    f = w(j) / ta(j);
    goto statement_55;
    statement_53:
    f = ta(j) / w(j);
    ta(j) = w(j);
    t = ta(jp1);
    ta(jp1) = tb(j);
    tb(j) = t;
    p(j) = tb(jp1);
    tb(jp1) = 0.e0;
    t = y(j);
    y(j) = y(jp1);
    y(jp1) = t;
    statement_55:
    tb(jp1) = tb(jp1) - f * p(j);
    ta(jp1) = ta(jp1) - f * tb(j);
    y(jp1) = y(jp1) - f * y(j);
    j = jp1;
  }
  if (ta(n) == 0.e0) {
    ta(n) = eps;
  }
  if (ta(l) == 0.e0) {
    ta(l) = eps;
  }
  y(n) = y(n) / ta(n);
  y(l) = (y(l) - y(n) * tb(l)) / ta(l);
  FEM_DO_SAFE(j, 2, l) {
    k = n - j;
    if (ta(k) == 0.e0) {
      ta(k) = eps;
    }
    y(k) = (y(k) - y(k + 1) * tb(k) - y(k + 2) * p(k)) / ta(k);
  }
  ay = fem::dabs(y(1));
  FEM_DO_SAFE(j, 2, n) {
    ay = fem::dmax1(ay, fem::dabs(y(j)));
  }
  FEM_DO_SAFE(j, 1, n) {
    y(j) = y(j) / ay;
  }
  ii++;
  if (ii <= 2) {
    goto statement_50;
  }
  id = ndim - 2;
  l = n - 2;
  FEM_DO_SAFE(j, 1, l) {
    id = id - j - 2;
    m = n - j;
    h = w(m - 1);
    if (h == 0.e0) {
      goto statement_68;
    }
    jp1 = j + 1;
    t = dot(c(id + 1), 1, y(m), 1, jp1) / (h * c(id + 1));
    kj = id;
    FEM_DO_SAFE(k, m, n) {
      kj++;
      y(k) += t * c(kj);
    }
    statement_68:;
  }
  xnorm = fem::dsqrt(dot(y, 1, y, 1, n));
  FEM_DO_SAFE(j, 1, n) {
    y(j) = y(j) / xnorm;
  }
  {
    write_loop wloop(cmn, iu, fem::unformatted);
    wloop, ev(ik);
    FEM_DO_SAFE(j, 1, n) {
      wloop, y(j);
    }
  }
  //C     WRITE(6,901) IK,EV(IK)
  ik++;
  if (ik <= n) {
    goto statement_1000;
  }
  cmn.io.endfile(iu);
  cmn.io.rewind(iu);
  //C
}

void
program_opendap_convert(
  int argc,
  char const* argv[])
{
  common cmn(argc, argv);
  common_read read(cmn);
  common_write write(cmn);
  fem::str<80> xyfl = fem::char0;
  fem::str<80> afl = fem::char0;
  fem::str<80> evcfl = fem::char0;
  int lmax = fem::int0;
  const int lmx = 40;
  int leny = fem::int0;
  int i = fem::int0;
  const int mxleny = fem::pow2((lmx + 1));
  const int mxmat = (mxleny * (mxleny + 1)) / 2;
  arr<double> ata(dimension(mxmat), fem::fill0);
  arr<float> wnorm(dimension(mxleny), fem::fill0);
  int nread = fem::int0;
  float xlon = fem::float0;
  float xlat = fem::float0;
  arr<float> y(dimension(mxleny), fem::fill0);
  const int mxwork = (lmx + 1) * 4;
  arr_1d<mxwork, float> wk1(fem::fill0);
  arr_1d<mxwork, float> wk2(fem::fill0);
  arr_1d<mxwork, float> wk3(fem::fill0);
  int k = fem::int0;
  int ind = fem::int0;
  int jj = fem::int0;
  int ii = fem::int0;
  arr<double> d1(dimension(mxleny), fem::fill0);
  arr<double> d2(dimension(mxleny), fem::fill0);
  arr<double> d3(dimension(mxleny), fem::fill0);
  arr<double> d4(dimension(mxleny), fem::fill0);
  arr<double> d5(dimension(mxleny), fem::fill0);
  arr<double> d6(dimension(mxleny), fem::fill0);
  arr<double> d7(dimension(mxleny), fem::fill0);
  arr<double> d8(dimension(mxleny), fem::fill0);
  arr<double> eigv(dimension(mxleny), fem::fill0);
  fem::str<80> ofl = fem::char0;
  float dum = fem::float0;
  double damp = fem::double0;
  float grdv = fem::float0;
  float xlon2 = fem::float0;
  float xlat2 = fem::float0;
  arr<float> atd(dimension(mxleny), fem::fill0);
  arr<double> evc(dimension(mxleny), fem::fill0);
  double sum = fem::double0;
  int j = fem::int0;
  double f1 = fem::double0;
  double w = fem::double0;
  arr<double> x(dimension(mxleny), fem::fill0);
  float dep1 = fem::float0;
  float dep2 = fem::float0;
  fem::str<80> sphfl = fem::char0;
  int id1 = fem::int0;
  int id2 = fem::int0;
  int np = fem::int0;
  arr<float> rawvec(dimension(mxleny), fem::fill0);
  float rcmb = fem::float0;
  float rmoho = fem::float0;
  float rearth = fem::float0;
  float dep = fem::float0;
  float r = fem::float0;
  const int maxp = 2891;
  int m = fem::int0;
  int mp = fem::int0;
  int ip = fem::int0;
  float xd = fem::float0;
  arr<float, 2> z(dimension(maxp + 1, 21), fem::fill0);
  const int maxd = 21;
  int ip1 = fem::int0;
  int ip2 = fem::int0;
  arr<float> d(dimension(maxp + 1), fem::fill0);
  float t = fem::float0;
  double damp2 = fem::double0;
  arr_1d<maxd, double> x2(fem::fill0);
  int idp1 = fem::int0;
  int idp2 = fem::int0;
  fem::str<120> strt = fem::char0;
  int lstrt = fem::int0;
  float spl = fem::float0;
  int ind1 = fem::int0;
  int l = fem::int0;
  float v = fem::float0;
  //C
  //C     construct matrices to expand data on a nxn degree grid
  //C     to spherical harmonics
  //C     Jeroen Ritsema
  //Cp    PK obtained from JR, 2013
  //C
  //C   ----------------------------------------------
  //C               invexpandxy.f
  //C   ----------------------------------------------
  //C
  //C   ----------------------------------------------
  //C                   sphexp.f
  //C   ----------------------------------------------
  //C
  //C     call chekcl('| :r:1:input x,y(,z) file (order= lon, lat!!!)'
  //C    1          //'|-a:r:1:output .a file'
  //C    1          //'|-evc:r:1:output .evc file'
  //C    1          //'|-lmax:r:1:maximum degree and order for expansion'
  //C    1          //'|')
  //C
  read(5, star), xyfl;
  read(5, star), afl;
  read(5, star), evcfl;
  read(5, star), lmax;
  //C
  //C      read(5,*) grdfl
  //C      read(5,*) ofl
  //C      read(5,*) afl
  //C      read(5,*) evcfl
  //C
  if (lmax > lmx) {
    FEM_STOP("lmx.gt.LMX");
  }
  leny = fem::pow2((lmax + 1));
  //C
  FEM_DO_SAFE(i, 1, mxmat) {
    ata(i) = 0.f;
  }
  //C
  cmn.io.open(19, xyfl)
    .status("old");
  cmn.io.open(21, afl)
    .form("unformatted")
    .status("unknown");
  cmn.io.open(22, evcfl)
    .form("unformatted")
    .status("unknown");
  //C
  //C     get normalisation
  normylm(lmax, wnorm);
  //C
  nread = 0;
  statement_10:
  try {
    read(19, star), xlon, xlat;
  }
  catch (fem::read_end const&) {
    goto statement_100;
  }
  nread++;
  if (fem::mod(nread, 1000) == 0) {
    write(6, "(i8,' points read')"), nread;
  }
  ylm(xlat, xlon, lmax, y, wk1, wk2, wk3);
  FEM_DO_SAFE(k, 1, leny) {
    y(k) = y(k) * wnorm(k);
  }
  //C
  ind = 0;
  FEM_DO_SAFE(jj, 1, leny) {
    FEM_DO_SAFE(ii, jj, leny) {
      ind++;
      ata(ind) += fem::dble(y(ii) * y(jj));
    }
  }
  //C
  {
    write_loop wloop(cmn, 21, fem::unformatted);
    wloop, xlon, xlat;
    FEM_DO_SAFE(k, 1, leny) {
      wloop, y(k);
    }
  }
  goto statement_10;
  //C
  statement_100:
  //C
  write(6, star), "decomposing......";
  write(22, fem::unformatted), lmax;
  ahouse2(leny, 22, ata, d1, d2, d3, d4, d5, d6, d7, d8, eigv);
  //C
  //C   ------------invexpandxy---------------
  //C
  //C       19 = inpm
  //C       21 = inpm.a
  //C       22 = inpm.evc
  //C
  //C      grdfl - inpm
  //C      ofl - inpm.raw
  //C      afl - inpm.a
  //C      evcfl -inpm.evc
  //C
  cmn.io.rewind(19);
  cmn.io.rewind(21);
  cmn.io.rewind(22);
  //C
  read(5, star), ofl;
  read(5, star), dum;
  damp = fem::dble(dum);
  //C
  read(22, fem::unformatted), lmax;
  leny = fem::pow2((lmax + 1));
  //C
  statement_12:
  try {
    read(19, star), xlon, xlat, grdv;
  }
  catch (fem::read_end const&) {
    goto statement_101;
  }
  {
    read_loop rloop(cmn, 21, fem::unformatted);
    rloop, xlon2, xlat2;
    FEM_DO_SAFE(k, 1, leny) {
      rloop, y(k);
    }
  }
  if (xlon != xlon2) {
    FEM_STOP("inconsistent .a and .xyz file");
  }
  if (xlat != xlat2) {
    FEM_STOP("inconsistent .a and .xyz file");
  }
  FEM_DO_SAFE(k, 1, leny) {
    atd(k) += grdv * y(k);
  }
  goto statement_12;
  //C
  statement_101:
  //C
  FEM_DO_SAFE(i, 1, leny) {
    {
      read_loop rloop(cmn, 22, fem::unformatted);
      rloop, eigv(i);
      FEM_DO_SAFE(k, 1, leny) {
        rloop, evc(k);
      }
    }
    //C
    if (eigv(i) > 1.e-7 * eigv(1)) {
      sum = 0.f;
      FEM_DO_SAFE(j, 1, leny) {
        sum += fem::dble(atd(j)) * evc(j);
      }
      //C
      f1 = 1.e0 / (eigv(i) + damp);
      w = sum * f1;
      FEM_DO_SAFE(j, 1, leny) {
        x(j) += w * evc(j);
      }
    }
  }
  //C
  //C     normaliseer harmonics
  normylm(lmax, wnorm);
  FEM_DO_SAFE(i, 1, leny) {
    x(i) = x(i) * fem::dble(wnorm(i));
  }
  //C
  cmn.io.open(24, ofl)
    .status("unknown");
  write(24, "(i3)"), lmax;
  //C-- Hendrik multiplied by 0.01 .... WHY?
  {
    write_loop wloop(cmn, 24, "(5e16.8)");
    FEM_DO_SAFE(i, 1, leny) {
      wloop, x(i) * .01f;
    }
  }
  //C--   write(24,'(5e16.8)') (x(i),i=1,leny)
  //C      close(24)
  //C
  //C   ------------sphexp---------------
  //C
  //C      rawfl - inpm.raw
  //C      dep1 - $dep1
  //C      dep2 - $dep2
  //C      sphfl - inpm.$iz.sph
  cmn.io.rewind(24);
  read(5, star), dep1;
  read(5, star), dep2;
  read(5, star), sphfl;
  //C
  if (dep2 < dep1) {
    FEM_STOP("STOP dep2 is less than dep1");
  }
  //C
  id1 = fem::nint(2891.f - dep2);
  id2 = fem::nint(2891.f - dep1);
  //C
  //C --  read raw file
  write(6, star), "model lmax =", lmax;
  np = fem::pow2((lmax + 1));
  if (lmax > lmx) {
    FEM_STOP("STOP >>>> lmax.gt.LMX");
  }
  write(6, star), "model np =", np;
  {
    read_loop rloop(cmn, 24, "(5e16.8)");
    FEM_DO_SAFE(i, 1, np) {
      rloop, rawvec(i);
    }
  }
  cmn.io.close(24);
  //C
  //C    Calculate the spline basis functions at a regular grid
  splhsetup(cmn);
  //C
  rcmb = 3480.f;
  rmoho = 6346.f;
  rearth = 6371.f;
  r = rearth - dep;
  m = maxp / 2;
  mp = maxp + 1;
  //C     xd=-1.+2.*(r-rcmb)/(rmoho-rcmb)
  FEM_DO_SAFE(ip, 1, 21) {
    FEM_DO_SAFE(j, -m, m) {
      xd = fem::ffloat(j) / fem::ffloat(m);
      z(j + m + 1, ip) = splh(cmn, ip - 1, xd);
      //C       write(20+ip,*) xd,z(j+m+1,ip)
    }
  }
  //C
  FEM_DO_SAFE(i, 1, mxmat) {
    ata(i) = 0.e0;
  }
  //C
  leny = maxd;
  //C
  FEM_DO_SAFE(ii, 1, mp) {
    ind = 0;
    FEM_DO_SAFE(ip1, 1, leny) {
      FEM_DO_SAFE(ip2, ip1, leny) {
        ind++;
        ata(ind) += fem::dble(z(ii, ip1) * z(ii, ip2));
      }
    }
  }
  //C
  cmn.io.open(101, "tmp.evc")
    .form("unformatted")
    .status("unknown");
  write(6, star), "decomposing......";
  ahouse(cmn, leny, 101, ata, d1, d2, d3, d4, d5, d6, d7, d8, eigv);
  cmn.io.close(101);
  //C
  cmn.io.open(101, "tmp.evc")
    .form("unformatted")
    .status("unknown");
  //C
  FEM_DO_SAFE(i, id1, id2) {
    d(i) = 1.f;
  }
  FEM_DO_SAFE(ip, 1, leny) {
    t = 0.f;
    FEM_DO_SAFE(j, 1, mp) {
      t += d(j) * z(j, ip);
    }
    atd(ip) = fem::dble(t);
  }
  //C
  FEM_DO_SAFE(i, 1, leny) {
    {
      read_loop rloop(cmn, 101, fem::unformatted);
      rloop, eigv(i);
      FEM_DO_SAFE(k, 1, leny) {
        rloop, evc(k);
      }
    }
    //C
    if (eigv(i) > 1.e-7 * eigv(1)) {
      sum = 0.f;
      FEM_DO_SAFE(j, 1, leny) {
        sum += fem::dble(atd(j)) * evc(j);
      }
      //C
      f1 = 1.e0 / (eigv(i) + damp2);
      w = sum * f1;
      FEM_DO_SAFE(j, 1, leny) {
        x2(j) += w * evc(j);
      }
    }
  }
  FEM_DO_SAFE(i, 1, leny) {
    write(6, star), "AA", i, x2(i);
  }
  //C
  ind = 1;
  cmn.io.open(26, sphfl)
    .status("unknown");
  idp1 = 4;
  idp2 = 24;
  write(6, star), "lmax,idp1,idp2=", lmax, idp1, idp2;
  wsphhead(cmn, lmax, idp1, idp2, strt, lstrt);
  write(26, star), strt(1, lstrt);
  FEM_DO_SAFE(j, idp1, idp2) {
    spl = x2(j - 3);
    write(6, star), "spl= ", spl;
    ind = 1;
    FEM_DO_SAFE(k, 0, lmax) {
      ind1 = ind + 2 * k;
      {
        write_loop wloop(cmn, 26, "(11e12.4)");
        FEM_DO_SAFE(l, ind, ind1) {
          wloop, spl * rawvec(l);
        }
      }
      ind = ind1 + 1;
    }
  }
  cmn.io.close(26);
  //C
  FEM_DO_SAFE(ii, 1, mp) {
    v = 0.f;
    FEM_DO_SAFE(ip, 1, leny) {
      v += z(ii, ip) * fem::sngl(x2(ip));
    }
    write(103, star), ii, v;
  }
  //C
}

} // namespace opendap_convert

int
main(
  int argc,
  char const* argv[])
{
  return fem::main_with_catch(
    argc, argv,
    opendap_convert::program_opendap_convert);
}
