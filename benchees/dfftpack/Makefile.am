if DOUBLE
if HAVE_F77
PRG=doit
endif
endif

INCLUDES = $(INCLBENCH)
EXTRA_DIST=test.f

FFTPACK=cfftb1.f cfftf1.f cffti1.f cosqb1.f cosqf1.f dcosqb.f dcosqf.f	\
dcosqi.f dcost.f dcosti.f dfftb.f dfftf.f dffti.f dsinqb.f dsinqf.f	\
dsinqi.f dsint.f dsinti.f dzfftb.f dzfftf.f dzffti.f ezfft1.f passb.f	\
passb2.f passb3.f passb4.f passb5.f passf.f passf2.f passf3.f passf4.f	\
passf5.f radb2.f radb3.f radb4.f radb5.f radbg.f radf2.f radf3.f	\
radf4.f radf5.f radfg.f rfftb1.f rfftf1.f rffti1.f sint1.f zfftb.f	\
zfftf.f zffti.f

doit_SOURCES=doit.c $(FFTPACK)
doit_LDADD=$(LIBBENCH) @FLIBS@

include ../Makefile.common

