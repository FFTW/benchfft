PKG=fft
TARFILE=$(PKG).tgz
PKG2=fft2d
TARFILE2=$(PKG2).tgz

if DOUBLE
if HAVE_F77
FORTRAN_PRG=doit_8gf doit_4gf doit_sgf doit_sg2df doit_4f2df doit_sg3df
FORTRAN_LIBS=lib8gf.a lib4gf.a libsgf.a libsg2df.a lib4f2df.a libsg3df.a
LOCAL_COPIES=fft4gf.f fft8gf.f fftsgf.f fftsg2df.f fft4f2df.f fftsg3df.f extract2 $(PKG2)
else
LOCAL_COPIES=extract2 $(PKG2)
endif
else
LOCAL_COPIES=extract2 $(PKG2)
endif

PRG=doit_sg doit_8g doit_4g doit_sg2d doit_4f2d doit_sg3d $(FORTRAN_PRG)
INCLUDES = $(INCLBENCH)
EXTRA_DIST = $(TARFILE) $(TARFILE2) doit.c doit-doc.h doit2d.c doit3d.c
noinst_LIBRARIES=$(FORTRAN_LIBS)

doit_sg_SOURCES=doit_sg.c
nodist_doit_sg_SOURCES=$(PKG)/fftsg.c
doit_sg_LDADD=$(LIBBENCH)
doit_sg.o: doit-doc.h doit.c

doit_8g_SOURCES=doit_8g.c
nodist_doit_8g_SOURCES=$(PKG)/fft8g.c
doit_8g_LDADD=$(LIBBENCH)
doit_8g.o: doit-doc.h doit.c

doit_4g_SOURCES=doit_4g.c
nodist_doit_4g_SOURCES=$(PKG)/fft4g.c
doit_4g_LDADD=$(LIBBENCH)
doit_4g.o: doit-doc.h doit.c

doit_sgf_SOURCES=doit_sgf.c
nodist_libsgf_a_SOURCES=fftsgf.f
doit_sgf_LDADD=libsgf.a $(LIBBENCH) @FLIBS@
doit_sgf.o: doit-doc.h doit.c

doit_8gf_SOURCES=doit_8gf.c
nodist_lib8gf_a_SOURCES=fft8gf.f
doit_8gf_LDADD=lib8gf.a $(LIBBENCH) @FLIBS@
doit_8gf.o: doit-doc.h doit.c

doit_4gf_SOURCES=doit_4gf.c
nodist_lib4gf_a_SOURCES=fft4gf.f
doit_4gf_LDADD=lib4gf.a $(LIBBENCH) @FLIBS@
doit_4gf.o: doit-doc.h doit.c

doit_sg2d_SOURCES=doit_sg2d.c
nodist_doit_sg2d_SOURCES=$(PKG)/fftsg.c $(PKG2)/fftsg2d.c
doit_sg2d_LDADD=$(LIBBENCH)
doit_sg2d.o: doit-doc.h doit2d.c

doit_sg2df_SOURCES=doit_sg2df.c
nodist_libsg2df_a_SOURCES=fftsgf.f fftsg2df.f
doit_sg2df_LDADD=libsg2df.a $(LIBBENCH) @FLIBS@
doit_sg2df.o: doit-doc.h doit2d.c

doit_4f2d_SOURCES=doit_4f2d.c
nodist_doit_4f2d_SOURCES=$(PKG2)/fft4f2d.c
doit_4f2d_LDADD=$(LIBBENCH)
doit_4f2d.o: doit-doc.h doit2d.c

doit_4f2df_SOURCES=doit_4f2df.c
nodist_lib4f2df_a_SOURCES=fft4f2df.f
doit_4f2df_LDADD=lib4f2df.a $(LIBBENCH) @FLIBS@
doit_4f2df.o: doit-doc.h doit2d.c

doit_sg3d_SOURCES=doit_sg3d.c
nodist_doit_sg3d_SOURCES=$(PKG)/fftsg.c $(PKG2)/fftsg3d.c
doit_sg3d_LDADD=$(LIBBENCH)
doit_sg3d.o: doit-doc.h doit3d.c

doit_sg3df_SOURCES=doit_sg3df.c
nodist_libsg3df_a_SOURCES=fftsgf.f fftsg3df.f
doit_sg3df_LDADD=libsg3df.a $(LIBBENCH) @FLIBS@
doit_sg3df.o: doit-doc.h doit3d.c

include ../Makefile.common

$(PKG)/fftsg.c $(PKG)/fft8g.c $(PKG)/fft4g.c: extract
$(PKG2)/fftsg2d.c $(PKG2)/fft4f2d.c $(PKG2)/fftsg3d.c: extract2

fftsgf.f: extract
	@CP@ -f $(PKG)/fftsg.f $@

fft8gf.f: extract
	@CP@ -f $(PKG)/fft8g.f $@

fft4gf.f: extract
	@CP@ -f $(PKG)/fft4g.f $@

fftsg2df.f: extract2
	@CP@ -f $(PKG2)/fftsg2d.f $@

fft4f2df.f: extract2
	@CP@ -f $(PKG2)/fft4f2d.f $@

fftsg3df.f: extract2
	@CP@ -f $(PKG2)/fftsg3d.f $@

# extract TARFILE2
extract2: $(TARFILE2)
	(@GZIP@ -dc ${srcdir}/$(TARFILE2) | @TAR@ xf -) && \
        echo timestamp >$@
