dnl Process this file with autoconf to produce a configure script.
AC_INIT(libbench/bench.h)
AC_PREREQ(2.50)
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE(benchfft, 3.0-alpha0)
AM_CONFIG_HEADER(config.h)

# configure options
AC_ARG_ENABLE(debug, [AC_HELP_STRING([--enable-debug],[compile benchfft with extra runtime checks for debugging])], ok=$enableval, ok=no)
if test "$ok" = "yes"; then
	AC_DEFINE(BENCHFFT_DEBUG,1,[Define to enable extra benchfft debugging code.])
fi

AC_ARG_ENABLE(single, [AC_HELP_STRING([--enable-single],[compile single-precision benchmark])], benchfft_single=$enableval, benchfft_single=no)
if test "$benchfft_single" = "yes"; then
	AC_DEFINE(BENCHFFT_SINGLE,1,[Define to benchmark in single precision.])
fi
AC_ARG_WITH(bothprec, [AC_HELP_STRING([--with-bothprec],[compile doits for both precisions])], benchfft_bothprec=$enableval, benchfft_bothprec=no)

AM_CONDITIONAL(SINGLE, test "$benchfft_single" = "yes" || test "$benchfft_bothprec" = "yes")
AM_CONDITIONAL(DOUBLE, test "$benchfft_single" = "no" || test "$benchfft_bothprec" = "yes")

# Checks for programs.
AC_PROG_CC
AC_PROG_CC_STDC
ACX_PROG_CC_EGCS
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_AWK

# check for c++, but allow the check to fail
AC_CHECK_PROGS(CXX, $CCC c++ g++ gcc CC cxx cc++ cl, nosuchprogram)
if test "$CXX" = nosuchprogram; then
    CXX=""
else
    AC_PROG_CXX
fi

AC_CHECK_PROG(GZIP, gzip, gzip)
AC_CHECK_PROG(TAR, tar, tar)
AC_CHECK_PROG(MKDIR, mkdir, mkdir)
AC_CHECK_PROG(TR, tr, tr)
AC_CHECK_PROG(CP, cp, cp)
AC_CHECK_PROG(PATCH, patch, patch)

# Until autoconf gets separate F90 support, we must require that *all*
# .f files be compiled F90 if we want anything to be compiled with F90.
AC_ARG_WITH(f90, [AC_HELP_STRING([--with-f90],[use F90 compiler for Fortran])],
	ok=$withval, ok=no)
if test "x$ok" = "xyes"; then
	AC_CHECK_PROGS(F90, f90 xlf90 pgf90 fort f95 xlf95 g95)
	test -n "$F90" && F77="$F90"
fi

AC_PROG_F77
if test -n "$F77"; then
	dnl replace with AC_F77_LIBRARY_LDFLAGS when autoconf is fixed. 
	AC_F77_LIBRARY_LDFLAGS
	AC_F77_DUMMY_MAIN([], bench_f77_fails=yes)
	if test "x$bench_f77_fails" = xyes; then
		AC_MSG_WARN([Cannot figure out how to link C code with Fortran libraries.])
		F77=""  # disable F77 benchmarks
		F90=""  # F90 too, since they use the same compiler
	else
		AC_F77_WRAPPERS
	fi
fi

test "x$F90FLAGS" = x && F90FLAGS="$FFLAGS"
test "x$F90LIBS" = x && F90LIBS="$FLIBS"
AC_SUBST(F90FLAGS)
AC_SUBST(F90LIBS)

AM_CONDITIONAL(HAVE_CXX, test -n "$CXX")
AM_CONDITIONAL(HAVE_F77, test -n "$F77")
AM_CONDITIONAL(HAVE_F90, test -n "$F90")

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS(malloc.h sys/time.h unistd.h string.h)
AC_CHECK_HEADERS(alloca.h stddef.h libintl.h)

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_HEADER_TIME

# Checks for library functions.
AC_FUNC_ALLOCA
AC_FUNC_MALLOC
AC_FUNC_VPRINTF
AC_CHECK_LIB(m, sin)
AC_CHECK_FUNCS(BSDgettimeofday)
AC_CHECK_FUNCS(gettimeofday)
AC_CHECK_FUNCS(drand48)
AC_CHECK_FUNCS(hypot)
AC_CHECK_FUNCS(sqrt)
AC_CHECK_FUNCS(strtod)
AC_CHECK_FUNCS(memset)
AC_CHECK_FUNCS(memalign)
AC_CHECK_DECLS(drand48)
AC_CHECK_DECLS(hypot,[],[],[#include <math.h>])

# Check for installed libraries
# check for pthreads library (required by MKL below)
AC_CHECK_LIB(pthread, pthread_attr_init)

# Check for FFTW 2
AC_CHECK_LIB(fftw, fftw_create_plan_specific,
	[LFFTW="-lfftw"; have_fftw2=yes])
AC_CHECK_LIB(rfftw, rfftw_create_plan_specific,
	[LFFTW="-lrfftw $LFFTW"; have_rfftw2=yes],[],$LFFTW)
AC_CHECK_HEADERS(fftw.h, [have_fftw_h=yes])
AC_CHECK_HEADERS(rfftw.h, [have_rfftw_h=yes])
AC_SUBST(LFFTW)
test -n "$have_fftw2" && test -n "$have_rfftw2" && test -n "$have_fftw_h" && test -n "$have_rfftw_h" && have_fftw2_all=yes
AM_CONDITIONAL(HAVE_FFTW2, test -n "$have_fftw2_all")
AC_CHECK_LIB(sfftw, fftw_create_plan_specific,
	[LSFFTW="-lsfftw"; have_sfftw2=yes])
AC_CHECK_LIB(srfftw, rfftw_create_plan_specific,
	[LSFFTW="-lsrfftw $LSFFTW"; have_srfftw2=yes],[],$LSFFTW)
AC_CHECK_HEADERS(sfftw.h)
AC_CHECK_HEADERS(srfftw.h)
AC_SUBST(LSFFTW)
test -n "$have_sfftw2" && test -n "$have_srfftw2" && have_sfftw2_all=yes
AM_CONDITIONAL(HAVE_SFFTW2, test -n "$have_sfftw2_all")

# Check for DEC Alpha DXML library
AC_CHECK_LIB(dxml, zfft_, [have_dxml=yes], [AC_CHECK_LIB(cxml, zfft_, [have_dxml=yes])])
AM_CONDITIONAL(HAVE_DXML, test -n "$have_dxml")
AC_CHECK_HEADERS(dxmldef.h cxmldef.h)

# Check for Sun Performance library:
AC_CHECK_LIB(fsu, __f90_stop)
SAVE_LIBS=$LIBS
AC_CHECK_LIB(sunperf, cfftf, [LIBS="-lsunperf $LIBS";  have_sunperf=yes])
AM_CONDITIONAL(HAVE_SUNPERF, test -n "$have_sunperf")
AC_CHECK_HEADERS(sunperf.h)
AC_CHECK_FUNCS(sunperf_version_)
LIBS=$SAVE_LIBS  # sunperf's symbols may interfere with fftpack.

# Check for SGIMATH library:
AC_CHECK_LIB(complib.sgimath, zfft1d, [have_sgimath=yes])
AM_CONDITIONAL(HAVE_SGIMATH, test -n "$have_sgimath")
AC_CHECK_HEADERS(fft.h)

# Check for ESSL library (some versions have only one precision):
if test "x$benchfft_single" = xyes; then
	essl_prefix=s
else
	essl_prefix=d
fi
AC_CHECK_LIB(essl, ${essl_prefix}cft, [have_essl=yes], [], $FLIBS)
AM_CONDITIONAL(HAVE_ESSL, test -n "$have_essl")

# Check for NAG library:
AC_F77_FUNC(c06frf)
if test "x$LIBNAG" = x; then
	AC_CHECK_LIB(nag, $c06frf, [have_nag=yes; LIBNAG=-lnag], [], $FLIBS)
else # user-specified $LIBNAG
	my_LIBS="$LIBS"
	LIBS="$LIBNAG $LIBS $FLIBS"
	AC_CHECK_FUNC($c06frf, [have_nag=yes])
	LIBS="$my_LIBS"
fi
AM_CONDITIONAL(HAVE_NAG, test -n "$have_nag")
AC_SUBST(LIBNAG)

# Check for IMSL library:
if test "x$benchfft_single" = xyes; then
	AC_F77_FUNC(f2tcf, imsl_f2tcf)
else
	AC_F77_FUNC(df2tcf, imsl_f2tcf)
fi
if test "x$LIBIMSL" = x; then
        AC_CHECK_LIB(imsl, $imsl_f2tcf, [have_imsl=yes; LIBIMSL=-limsl], [], $FLIBS)
else # user-specified $LIBIMSL
        my_LIBS="$LIBS"
        LIBS="$LIBIMSL $LIBS $FLIBS"
        AC_CHECK_FUNC($imsl_f2tcf, [have_imsl=yes])
        LIBS="$my_LIBS"
fi
AM_CONDITIONAL(HAVE_IMSL, test -n "$have_imsl")
AC_SUBST(LIBIMSL)

# Check for Apple's vDSP/vecLib library:
xsave_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="-I/System/Library/Frameworks/vecLib.framework/Headers"
test -z "$VECLIB" && VECLIB=/System/Library/Frameworks/vecLib.framework/vecLib
AC_CHECK_HEADER(vDSP.h,
	[AC_CHECK_LIB(m, fft_zop, [have_vdsp=yes],
		[CPPFLAGS="$xsave_CPPFLAGS"], $VECLIB)],
	[CPPFLAGS="$xsave_CPPFLAGS"])
AC_SUBST(VECLIB)
AM_CONDITIONAL(HAVE_VDSP, test -n "$have_vdsp")
if test -n "$have_vdsp"; then
	xsave_LIBS="$LIBS"
	LIBS="$LIBS $VECLIB"
	AC_CHECK_FUNCS(vec_malloc, vec_free)
	LIBS="$xsave_LIBS"
fi

# Check for Intel's MKL library:
AC_CHECK_LIB(guide, omp_in_parallel)
AC_CHECK_LIB(mkl_def, zfft1d_, [have_mkl_def=yes])
AM_CONDITIONAL(HAVE_MKL_DEF, test -n "$have_mkl_def")
AC_CHECK_LIB(mkl_p3, zfft1d_, [have_mkl_p3=yes])
AM_CONDITIONAL(HAVE_MKL_P3, test -n "$have_mkl_p3")
AC_CHECK_LIB(mkl_p4, zfft1d_, [have_mkl_p4=yes])
AM_CONDITIONAL(HAVE_MKL_P4, test -n "$have_mkl_p4")
AC_CHECK_LIB(mkl_vml_itp, zfft1d_, [have_mkl_vml_itp=yes])
AM_CONDITIONAL(HAVE_MKL_VML_ITP, test -n "$have_mkl_vml_itp")
AC_CHECK_HEADERS(mkl_fft.h)
AC_CHECK_HEADERS(mkl_fftc_ln.h)
AC_CHECK_HEADERS(mkl_blas.h)
AC_CHECK_LIB(mkl_def, MKLGetVersionString, [AC_DEFINE(HAVE_MKLGETVERSIONSTRING, 1, "")])

# check for Numerical Recipes source in benchees/nr
nr_dir=$srcdir/benchees/nr
test -r $nr_dir/four1.c && test -r $nr_dir/fourn.c && test -r $nr_dir/realft.c && test -r $nr_dir/rlft3.c && have_nrc=yes
test -r $nr_dir/four1.f && test -r $nr_dir/fourn.f && test -r $nr_dir/realft.f && test -r $nr_dir/rlft3.f && have_nrf=yes

dnl Check for the GNU Scientific Library
SAVE_LIBS=$LIBS
AC_CHECK_LIB(gslcblas, cblas_dsdot) dnl why on the hell does it use cblas?
AC_CHECK_LIB(gsl, gsl_fft_complex_forward, [have_gsl=yes])
AC_CHECK_HEADERS([gsl/gsl_fft.h gsl/gsl_fft_complex.h gsl/gsl_fft_real.h])
AC_CHECK_HEADERS([gsl/gsl_fft_halfcomplex.h gsl/gsl_version.h])
AM_CONDITIONAL(HAVE_GSL, test -n "$have_gsl")
LIBS=$SAVE_LIBS

dnl Check for SCSL library (SGI/Cray Scientific Library):
AC_CHECK_LIB(scs, zzfft)

dnl Check for Cray SCILIB library:
AC_CHECK_LIB(sci, CCFFT)

AC_MSG_CHECKING([for isnan])
AC_TRY_LINK([#include <math.h>
], if (!isnan(3.14159)) isnan(2.7183);, ok=yes, ok=no)
if test "$ok" = "yes"; then
	AC_DEFINE(HAVE_ISNAN,1,[Define if the isnan() function/macro is available.])
fi
AC_MSG_RESULT(${ok})


dnl check for optimization options
dnl This macro sets ANSI mode for certain compilers, and must
dnl thus come before the subsequenct checks
ACX_PROG_CC_MAXOPT(benchfft)

dnl override CFLAGS selection when debugging
if test "${enable_debug}" = "yes"; then
	CFLAGS="-g"
fi

AC_DEFINE_UNQUOTED(BENCH_CC, "$CC $CFLAGS", [C compiler name and flags])
if test -n "$CXX"; then
	AC_DEFINE_UNQUOTED(BENCH_CXX, "$CXX $CXXFLAGS",
		[C++ compiler name and flags])
fi
if test -n "$F77"; then
	AC_DEFINE_UNQUOTED(BENCH_F77, "$F77 $FFLAGS",
		[Fortran 77 compiler name and flags])
fi
if test -n "$F90"; then
	AC_DEFINE_UNQUOTED(BENCH_F90, "$F90 $F90FLAGS",
		[Fortran 90 compiler name and flags])
fi

dnl TODO
ACX_GCC_ALIGNS_STACK()

AC_OUTPUT([
Makefile
libbench/Makefile
benchees/Makefile
scripts/Makefile

benchees/arprec/Makefile
benchees/athfft/Makefile
benchees/bloodworth/Makefile
benchees/burrus/Makefile
benchees/cross/Makefile
benchees/cwplib/Makefile
benchees/dfftpack/Makefile
benchees/djbfft-0.76/Makefile
benchees/dsp/Makefile
benchees/dxml/Makefile
benchees/emayer/Makefile
benchees/essl/Makefile
benchees/fftpack/Makefile
benchees/fftw2/Makefile
benchees/fxt/Makefile
benchees/glassman/Makefile
benchees/goedecker/Makefile
benchees/gpfa/Makefile
benchees/green-ffts-2.0/Makefile
benchees/gsl/Makefile
benchees/harm/Makefile
benchees/imsl/Makefile
benchees/intel-mkl/Makefile
benchees/jmfft/Makefile
benchees/krukar/Makefile
benchees/mfft/Makefile
benchees/mixfft/Makefile
benchees/monnier/Makefile
benchees/mpfun77/Makefile
benchees/mpfun90/Makefile
benchees/nag/Makefile
benchees/napack/Makefile
benchees/nr/Makefile
benchees/ooura/Makefile
benchees/qft/Makefile
benchees/ransom/Makefile
benchees/rmayer/Makefile
benchees/sciport/Makefile
benchees/sgimath/Makefile
benchees/singleton/Makefile
benchees/sorensen/Makefile
benchees/spiral-fft/Makefile
benchees/statlib/Makefile
benchees/sunperf/Makefile
benchees/temperton/Makefile
benchees/teneyck/Makefile
benchees/valkenburg/Makefile
benchees/vdsp/Makefile

])

cat <<EOF
*******************************************************************************
		    Optional FFT libraries Found:

                                       DXML/CXML (Alpha only): ${have_dxml-no}
                         Sun Performance Library (SPARC only): ${have_sunperf-no}
                                      sgimath (SGI/MIPS only): ${have_sgimath-no}
                                          ESSL (IBM/AIX only): ${have_essl-no}
                                  Apple VDSP (Macintosh only): ${have_vdsp-no}

                               NAG (proprietary, www.nag.com): ${have_nag-no}
                              IMSL (proprietary, www.vni.com): ${have_imsl-no}
             Numerical Recipes in C (proprietary, www.nr.com): ${have_nrc-no}
       Numerical Recipes in Fortran (proprietary, www.nr.com): ${have_nrf-no}

Intel MKL (gratis, www.intel.com/software/products/mkl/mkl52): ${have_mkl-no}

                                FFTW 2.x (free, www.fftw.org): ${have_fftw2_all-no}
               FFTW 2.x single-precision (free, www.fftw.org): ${have_sfftw2_all-no}
  GNU Scientific Library (GSL) (free, sources.redhat.com/gsl): ${have_gsl-no}
EOF
