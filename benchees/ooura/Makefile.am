if HAVE_F77
FORTRAN_PRG=doit_8gf doit_4gf
LOCAL_COPIES=fft4gf.f fft8gf.f
else
FORTRAN_PRG=
endif

PRG=doit_sg doit_8g doit_4g $(FORTRAN_PRG)
INCLUDES = $(INCLBENCH)
PKG=fft
TARFILE=$(PKG).tgz
EXTRA_DIST = $(TARFILE) doit.c

doit_sg_SOURCES=doit_sg.c
nodist_doit_sg_SOURCES=$(PKG)/fftsg.c
doit_sg_LDADD=$(LIBBENCH)
doit_sg.o: doit.c

doit_8g_SOURCES=doit_8g.c
nodist_doit_8g_SOURCES=$(PKG)/fft8g.c
doit_8g_LDADD=$(LIBBENCH)
doit_8g.o: doit.c

doit_4g_SOURCES=doit_4g.c
nodist_doit_4g_SOURCES=$(PKG)/fft4g.c
doit_4g_LDADD=$(LIBBENCH)
doit_4g.o: doit.c

doit_8gf_SOURCES=doit_8gf.c
nodist_doit_8gf_SOURCES=fft8gf.f
doit_8gf_LDADD=$(LIBBENCH) @FLIBS@
doit_8gf.o: doit.c

doit_4gf_SOURCES=doit_4gf.c
nodist_doit_4gf_SOURCES=fft4gf.f
doit_4gf_LDADD=$(LIBBENCH) @FLIBS@
doit_4gf.o: doit.c

include ../Makefile.common

$(PKG)/fftsg.c $(PKG)/fft8g.c $(PKG)/fft4g.c: extract

fft8gf.f: extract
	@CP@ -f $(PKG)/fft8g.f $@

fft4gf.f: extract
	@CP@ -f $(PKG)/fft4g.f $@
