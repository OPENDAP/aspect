#include <fem.hpp> // Fortran EMulation library of fable module

namespace mkexpmatxy {

using namespace fem::major_types;

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

using fem::common;

void
program_mkexpmatxy(
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
  //C
  //C     construct matrices to expand data on a nxn degree grid
  //C     to spherical harmonics
  //C     Jeroen Ritsema
  //Cp    PK obtained from JR, 2013
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
  // FIXME normylm(lmax, wnorm);
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
  // FIXME ylm(xlat, xlon, lmax, y, wk1, wk2, wk3);
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
  // FIXME(leny, 22, ata, d1, d2, d3, d4, d5, d6, d7, d8, eigv);
  //C
}

} // namespace mkexpmatxy

int
main(
  int argc,
  char const* argv[])
{
  return fem::main_with_catch(
    argc, argv,
    mkexpmatxy::program_mkexpmatxy);
}
