if SINGLE
if HAVE_F77
PRG=doit
noinst_LIBRARIES=libmfft.a
endif
endif

INCLUDES = $(INCLBENCH)

PKG = mfft
TARFILE = $(PKG).tar.gz
EXTRA_DIST = $(TARFILE)

doit_SOURCES=doit.c
doit_LDADD=libmfft.a $(LIBBENCH) @FLIBS@

MFFT_SRC = $(PKG)/c2fft.f $(PKG)/c3fft.f $(PKG)/c4fft.f		\
$(PKG)/mffta4.f $(PKG)/mffta5.f $(PKG)/mffta6.f $(PKG)/mffta7.f	\
$(PKG)/mffta8.f $(PKG)/mffta9.f $(PKG)/mfftb4.f $(PKG)/mfftb5.f	\
$(PKG)/mfftb6.f $(PKG)/mfftb7.f $(PKG)/mfftb8.f $(PKG)/mfftb9.f	\
$(PKG)/mfftc4.f $(PKG)/mfftc5.f $(PKG)/mfftc6.f $(PKG)/mfftc7.f	\
$(PKG)/mfftc8.f $(PKG)/mfftc9.f $(PKG)/mfftdm.f $(PKG)/mfftds.f	\
$(PKG)/mfftdv.f $(PKG)/mfftim.f $(PKG)/mfftis.f $(PKG)/mfftiv.f	\
$(PKG)/mfftom.f $(PKG)/mfftov.f $(PKG)/mfftp.f $(PKG)/mfftp1.f	\
$(PKG)/mfftp2.f $(PKG)/mfftp3.f $(PKG)/mfftp4.f $(PKG)/mfftrd.f	\
$(PKG)/mfftri.f $(PKG)/mfftrp.f $(PKG)/mfftz0.f $(PKG)/r2fft.f	\
$(PKG)/r3fft.f $(PKG)/r4fft.f

libmfft_a_SOURCES = $(MFFT_SRC)

$(MFFT_SRC): extract

include ../Makefile.common
