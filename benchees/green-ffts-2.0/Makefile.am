noinst_PROGRAMS=doit
INCLUDES = -I${top_srcdir}/libbench -Ddouble_t=double
PKG=green-ffts-2.0
TARFILE=$(PKG).tar.gz
EXTRA_DIST = $(TARFILE) fp.h

doit_SOURCES=doit.c
doit_LDADD=dxpose.o fft2d.o fftext.o fftlib.o matlib.o	\
${top_builddir}/libbench/libbench.a

doit.o: $(PKG)

$(PKG): $(TARFILE)
	gzip -dc ${srcdir}/$(TARFILE) | tar xf -

dxpose.c dxpose.h fft2d.c fft2d.h fftext.c fftext.h fftlib.c fftlib.h matlib.c matlib.h: $(PKG)
	tr '\r' '\n' <$(PKG)/fft\ libraries/$@ >$@

# explicit dependencies because source file need to be created
dxpose.o: dxpose.c dxpose.h
fft2d.o: fft2d.c fp.h fftlib.h fftext.h matlib.h dxpose.h fft2d.h
fftext.o: fftext.c fftlib.h matlib.h fftext.h
fftlib.o: fftlib.c fftlib.h
matlib.o: matlib.c matlib.h
doit.o: fftext.h