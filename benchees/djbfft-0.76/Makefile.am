noinst_PROGRAMS=doit
INCLUDES = -I${top_srcdir}/libbench -I./djbfft-0.76
EXTRA_DIST = djbfft-0.76.tar.gz

doit_SOURCES=doit.c
doit_LDADD=djbfft-0.76/djbfft.a ${top_builddir}/libbench/libbench.a

doit.o: djbfft-0.76

djbfft-0.76/djbfft.a: djbfft-0.76
	(cd djbfft-0.76; make)

djbfft-0.76: djbfft-0.76.tar.gz
	gzip -dc ${srcdir}/djbfft-0.76.tar.gz | tar xf -

djbfft-0.76/fftc4.h djbfft-0.76/fftc8.h djbfft-0.76/fftfreq.h: djbfft-0.76