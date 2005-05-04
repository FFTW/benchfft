if DOUBLE
PRG=doit
endif

PKG = scimark2c
TARFILE = $(PKG).tar.gz
EXTRA_DIST = $(TARFILE)

INCLUDES = $(INCLBENCH) -I$(PKG)

FFT = $(PKG)/FFT.c $(PKG)/FFT.h

doit_SOURCES=doit.c 
nodist_doit_SOURCES=$(FFT)
doit_LDADD=$(LIBBENCH)
doit.o: extract

$(FFT): extract

include ../Makefile.common
