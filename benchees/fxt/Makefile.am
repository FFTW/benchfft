if HAVE_CXX
PRG=doit_dif doit_dit doit_split doit_fht doit_4step
else
FORTRAN_PRG=
endif

INCLUDES = $(INCLBENCH)
PKG=fxt
TARFILE=$(PKG)-2001.06.16.tar.gz
EXTRA_DIST = $(TARFILE) doit.cc

CPPFLAGS = @CPPFLAGS@ -I$(PKG)/include -I$(PKG)/auxil -I$(PKG)/permute -I$(PKG)/auxbit

doit_dit_SOURCES = doit_dit.cc
doit_dit_LDADD = libfxt.a $(LIBBENCH)
doit_dit.cc: doit.cc

doit_dif_SOURCES = doit_dif.cc
doit_dif_LDADD = libfxt.a $(LIBBENCH)
doit_dif.cc: doit.cc

doit_split_SOURCES = doit_split.cc
doit_split_LDADD = libfxt.a $(LIBBENCH)
doit_split.cc: doit.cc

doit_fht_SOURCES = doit_fht.cc
doit_fht_LDADD = libfxt.a $(LIBBENCH)
doit_fht.cc: doit.cc

doit_4step_SOURCES = doit_4step.cc
doit_4step_LDADD = libfxt.a $(LIBBENCH)
doit_4step.cc: doit.cc

noinst_LIBRARIES = libfxt.a

libfxt_a_SOURCES = $(PKG)/fft/cfftdif4.cc $(PKG)/fft/cfftdit4.cc	\
$(PKG)/fft/cfftsplitradix.cc $(PKG)/fft/cfftwrap.cc			\
$(PKG)/fft/fft8difcore.cc $(PKG)/fft/fft8ditcore.cc $(PKG)/fft/fft9.cc	\
$(PKG)/fft/fftcnvl.cc $(PKG)/fft/fftcnvla.cc $(PKG)/fft/fftcocnvl.cc	\
$(PKG)/fft/fftcocorr.cc $(PKG)/fft/fftcorr.cc $(PKG)/fft/fftdif4.cc	\
$(PKG)/fft/fftdit4.cc $(PKG)/fft/fftspect.cc				\
$(PKG)/fft/fftsplitradix.cc $(PKG)/fft/fouriershift.cc			\
$(PKG)/fft/skipfft.cc $(PKG)/matrix/columnffts.cc			\
$(PKG)/matrix/ctranspos2.cc $(PKG)/matrix/ldn2rc.h			\
$(PKG)/matrix/matrixfft.cc $(PKG)/matrix/matrixfftcnvl.cc		\
$(PKG)/matrix/matrixfftcnvla.cc $(PKG)/matrix/matrixfftcocnvla.cc	\
$(PKG)/matrix/rowcnvls.cc $(PKG)/matrix/rowffts.cc			\
$(PKG)/matrix/transpos2.cc $(PKG)/fht/cfht0.cc $(PKG)/fht/cfhtcnvl.cc	\
$(PKG)/fht/cfhtcnvla.cc $(PKG)/fht/cfhtdif.cc $(PKG)/fht/cfhtdit.cc	\
$(PKG)/fht/fht0.cc $(PKG)/fht/fhtcfft.cc $(PKG)/fht/fhtcnvl.cc		\
$(PKG)/fht/fhtcnvla.cc $(PKG)/fht/fhtcorr.cc $(PKG)/fht/fhtdif.cc	\
$(PKG)/fht/fhtdit.cc $(PKG)/fht/fhtfft.cc $(PKG)/fht/fhtmayer.txt	\
$(PKG)/fht/fhtnegacnvla.cc $(PKG)/fht/fhtspect.cc			\
$(PKG)/fht/fhttable.cc $(PKG)/fht/hartleyshift.cc			\
$(PKG)/fht/shortfhtdifcore.h $(PKG)/fht/shortfhtditcore.h		\
$(PKG)/fht/skipfht.cc $(PKG)/fht/twodimfht.cc				\
$(PKG)/fht/twodimfhtcnvl.cc $(PKG)/realfft/realfftbyfht.cc		\
$(PKG)/realfft/realffteasyord.cc $(PKG)/realfft/realfftsplitradix.cc	\
$(PKG)/realfft/realfftwrap.cc $(PKG)/realfft/skiprealfft.cc		\
$(PKG)/auxil/auxdouble.cc $(PKG)/auxil/auxprint.cc $(PKG)/auxil/gcd.cc	\
$(PKG)/auxil/jjassert.cc $(PKG)/auxil/resample.cc			\
$(PKG)/auxil/sincos.cc $(PKG)/auxil/trigtable.cc			\
$(PKG)/auxil/version.cc $(PKG)/auxil/workspace.cc

include ../Makefile.common
