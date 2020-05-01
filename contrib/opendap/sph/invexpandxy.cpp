#include <fem.hpp> // Fortran EMulation library of fable module

namespace invexpandxy {

using namespace fem::major_types;

void
normylm(...)
{
  throw std::runtime_error(
    "Missing function implementation: normylm");
}

using fem::common;

void
program_invexpandxy(
  int argc,
  char const* argv[])
{
  common cmn(argc, argv);
  common_read read(cmn);
  common_write write(cmn);
  fem::str<80> grdfl = fem::char0;
  fem::str<80> ofl = fem::char0;
  fem::str<80> afl = fem::char0;
  fem::str<80> evcfl = fem::char0;
  float dum = fem::float0;
  double damp = fem::double0;
  int lmax = fem::int0;
  int leny = fem::int0;
  float xlon = fem::float0;
  float xlat = fem::float0;
  float grdv = fem::float0;
  float xlon2 = fem::float0;
  float xlat2 = fem::float0;
  const int mxl = 40;
  const int mxleny = fem::pow2((mxl + 1));
  arr<float> y(dimension(mxleny), fem::fill0);
  int k = fem::int0;
  arr<float> atd(dimension(mxleny), fem::fill0);
  int i = fem::int0;
  arr<double> eigv(dimension(mxleny), fem::fill0);
  arr<double> evc(dimension(mxleny), fem::fill0);
  double sum = fem::double0;
  int j = fem::int0;
  double f1 = fem::double0;
  double w = fem::double0;
  arr<double> x(dimension(mxleny), fem::fill0);
  arr<float> wnorm(dimension(mxleny), fem::fill0);
  //C
  //C     least squares expansion in spherical harmonics of data
  //C     on a grid
  //Cp    PK obtained from JR, 2013
  //Cp    removed multiplication factor of 0.01 in writing out the .xyz file
  //Cp    Paula Koelemeijer, February 2015
  //C
  //C     call chekcl('|  :r:2:Input xyz file, output .raw file'
  //C    1          //'|-a:r:1:input .a file'
  //C    1          //'|-evc:r:1:evc file'
  //C    1          //'|-damp:o:1:damping'
  //C    1          //'|')
  //C
  read(5, star), grdfl;
  read(5, star), ofl;
  read(5, star), afl;
  read(5, star), evcfl;
  read(5, star), dum;
  damp = fem::dble(dum);
  //C
  cmn.io.open(21, afl)
    .form("unformatted")
    .status("old");
  cmn.io.open(22, evcfl)
    .form("unformatted")
    .status("old");
  cmn.io.open(23, grdfl)
    .status("old");
  //C
  read(22, fem::unformatted), lmax;
  leny = fem::pow2((lmax + 1));
  //C
  statement_12:
  try {
    read(23, star), xlon, xlat, grdv;
  }
  catch (fem::read_end const&) {
    goto statement_100;
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
  statement_100:
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
  cmn.io.close(24);
  //C
}

} // namespace invexpandxy

int
main(
  int argc,
  char const* argv[])
{
  return fem::main_with_catch(
    argc, argv,
    invexpandxy::program_invexpandxy);
}
