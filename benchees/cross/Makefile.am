PRG=doit

PKG = fft
TARFILE = $(PKG).tar.gz
EXTRA_DIST = $(TARFILE)

AM_CPPFLAGS = $(INCLBENCH) -I$(PKG)

FFT = $(PKG)/fourierd.c $(PKG)/fourierf.c $(PKG)/fftmisc.c

doit_SOURCES=doit.c
nodist_doit_SOURCES=$(FFT)
doit_LDADD=$(LIBBENCH)
doit.o: extract

$(FFT): extract

include ../Makefile.common
