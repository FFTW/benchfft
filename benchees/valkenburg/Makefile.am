if DOUBLE
PRG=doit
endif

PKG = fft2
TARFILE = $(PKG).tar.gz
EXTRA_DIST = $(TARFILE)

INCLUDES = $(INCLBENCH) -I$(PKG)

FFT = $(PKG)/complex/ft.c $(PKG)/complex/w.c $(PKG)/complex/fourier.c

doit_SOURCES=doit.c
nodist_doit_SOURCES=$(FFT)
doit_LDADD=$(LIBBENCH) @FLIBS@

doit.c $(FFT): extract

include ../Makefile.common
