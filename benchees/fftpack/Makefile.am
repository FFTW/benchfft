if SINGLE
if HAVE_F77
PRG=doit
noinst_LIBRARIES=libfftpack.a
endif
endif

INCLUDES = $(INCLBENCH)
EXTRA_DIST=test.f

libfftpack_a_SOURCES=cfftb.f cfftb1.f cfftf.f cfftf1.f cffti.f		\
cffti1.f cosqb.f cosqb1.f cosqf.f cosqf1.f cosqi.f cost.f costi.f	\
ezfft1.f ezfftb.f ezfftf.f ezffti.f passb.f passb2.f passb3.f passb4.f	\
passb5.f passf.f passf2.f passf3.f passf4.f passf5.f radb2.f radb3.f	\
radb4.f radb5.f radbg.f radf2.f radf3.f radf4.f radf5.f radfg.f		\
rfftb.f rfftb1.f rfftf.f rfftf1.f rffti.f rffti1.f sinqb.f sinqf.f	\
sinqi.f sint.f sint1.f sinti.f

doit_SOURCES=doit.c
doit_LDADD=libfftpack.a $(LIBBENCH) @FLIBS@

include ../Makefile.common
