if DOUBLE
if HAVE_FC
PRG=doit
noinst_LIBRARIES=libfft.a
endif
endif

INCLUDES = $(INCLBENCH)
EXTRA_DIST = sgfft.tar.gz

libfft_a_SOURCES=basis_defs.f90 sg_ctrig.f90 sg_fft.f90 sg_fftpad.f90	\
sg_fftpx.f90 sg_fftrisc.f90 sg_fftx.f90 sg_ffty.f90 sg_fftz.f90

doit_SOURCES=doit.c
doit_LDADD=libfft.a $(LIBBENCH) @FCLIBS@

# make sure basis_defs is compiled first

basis_defs.o: basis_defs.f90
	$(FC) -c $(FCFLAGS) $(FCFLAGS_f90) ${srcdir}/basis_defs.f90

sg_ctrig.o: basis_defs.o sg_ctrig.f90
	$(FC) -c $(FCFLAGS) $(FCFLAGS_f90) ${srcdir}/sg_ctrig.f90

sg_fft.o: basis_defs.o sg_fft.f90
	$(FC) -c $(FCFLAGS) $(FCFLAGS_f90) ${srcdir}/sg_fft.f90

sg_fftpad.o: basis_defs.o sg_fftpad.f90
	$(FC) -c $(FCFLAGS) $(FCFLAGS_f90) ${srcdir}/sg_fftpad.f90

sg_fftpx.o: basis_defs.o sg_fftpx.f90
	$(FC) -c $(FCFLAGS) $(FCFLAGS_f90) ${srcdir}/sg_fftpx.f90

sg_fftrisc.o: basis_defs.o sg_fftrisc.f90
	$(FC) -c $(FCFLAGS) $(FCFLAGS_f90) ${srcdir}/sg_fftrisc.f90

sg_fftx.o: basis_defs.o sg_fftx.f90
	$(FC) -c $(FCFLAGS) $(FCFLAGS_f90) ${srcdir}/sg_fftx.f90

sg_ffty.o: basis_defs.o sg_ffty.f90
	$(FC) -c $(FCFLAGS) $(FCFLAGS_f90) ${srcdir}/sg_ffty.f90

sg_fftz.o: basis_defs.o sg_fftz.f90
	$(FC) -c $(FCFLAGS) $(FCFLAGS_f90) ${srcdir}/sg_fftz.f90

include ../Makefile.common
