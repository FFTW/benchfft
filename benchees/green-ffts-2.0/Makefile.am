PRG=doit
INCLUDES = $(INCLBENCH) -Ddouble_t=double
PKG=green-ffts-2.0
TARFILE=$(PKG).tar.gz
EXTRA_DIST = $(TARFILE) fp.h
LOCAL_COPIES=dxpose.c dxpose.h fft2d.c fft2d.h fftext.c fftext.h	\
fftlib.c fftlib.h matlib.c matlib.h

doit_SOURCES=doit.c
doit_LDADD=dxpose.o fft2d.o fftext.o fftlib.o matlib.o $(LIBBENCH)

include ../Makefile.common

doit.o: extract

# remove Macintosh newlines and C++-style  comments
$(LOCAL_COPIES): extract
	(echo "#include <config.h>"; tr '\r' '\n' <$(PKG)/fft\ libraries/$@; echo ""; echo "") | sed -e 's/\/\/.*$$//g' >$@

# explicit dependencies because source file need to be created
dxpose.o: dxpose.c dxpose.h
fft2d.o: fft2d.c fp.h fftlib.h fftext.h matlib.h dxpose.h fft2d.h
fftext.o: fftext.c fftlib.h matlib.h fftext.h
fftlib.o: fftlib.c fftlib.h
matlib.o: matlib.c matlib.h
doit.o: fftext.h fft2d.h