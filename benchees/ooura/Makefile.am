if HAVE_F77
FORTRAN_PRG=doit_8gf doit_4gf
else
FORTRAN_PRG=
endif

noinst_PROGRAMS=doit_sg doit_8g doit_4g $(FORTRAN_PRG)
INCLUDES = -I${top_srcdir}/libbench
PKG=fft
TARFILE=$(PKG).tgz
EXTRA_DIST = $(TARFILE)

doit_sg_SOURCES=
doit_sg_LDADD=doit_sg.o fftsg.o ${top_builddir}/libbench/libbench.a
EXTRA_doit_sg_SOURCES=doit.c

doit_8g_SOURCES=
doit_8g_LDADD=doit_8g.o fft8g.o ${top_builddir}/libbench/libbench.a

doit_4g_SOURCES=
doit_4g_LDADD=doit_4g.o fft4g.o ${top_builddir}/libbench/libbench.a

doit_8gf_SOURCES=
doit_8gf_LDADD=doit_8gf.o fft8gf.o ${top_builddir}/libbench/libbench.a
EXTRA_doit_8gf_SOURCES=dummy.f  # force automake to include FORTRAN support

doit_4gf_SOURCES=
doit_4gf_LDADD=doit_4gf.o fft4gf.o ${top_builddir}/libbench/libbench.a

extract: $(TARFILE)
	gzip -dc ${srcdir}/$(TARFILE) | tar xf -
	@echo timestamp >$@

doit_sg.o: doit.c $(CONFIG_HEADER)
	$(COMPILE) -c -o $@ -DNAME=\"ooura-sg\" ${srcdir}/doit.c

fftsg.o: extract $(CONFIG_HEADER)
	$(COMPILE) -c -o $@ $(PKG)/fftsg.c

doit_8g.o: doit.c $(CONFIG_HEADER)
	$(COMPILE) -c -o $@ -DNAME=\"ooura-8g\" ${srcdir}/doit.c

fft8g.o: extract
	$(COMPILE) -c -o $@ $(PKG)/fft8g.c 

doit_4g.o: doit.c $(CONFIG_HEADER)
	$(COMPILE) -c -o $@ -DNAME=\"ooura-4g\" ${srcdir}/doit.c

fft4g.o: extract
	$(COMPILE) -c -o $@ $(PKG)/fft4g.c

doit_8gf.o: doit.c $(CONFIG_HEADER)
	$(COMPILE) -c -o $@ -DNAME=\"ooura-8gf\" -DFORTRAN ${srcdir}/doit.c

fft8gf.o: extract
	$(F77COMPILE) -c -o $@ $(PKG)/fft8g.f

doit_4gf.o: doit.c $(CONFIG_HEADER)
	$(COMPILE) -c -o $@ -DNAME=\"ooura-4gf\" -DFORTRAN ${srcdir}/doit.c

fft4gf.o: extract
	$(F77COMPILE) -c -o $@ $(PKG)/fft4g.f
