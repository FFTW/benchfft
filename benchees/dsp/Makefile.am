if HAVE_F77
PRG=doit11 doit12 doit12-fast doit12-ffa
noinst_LIBRARIES=libdsp.a
endif

PKG=dsp
TARFILE=$(PKG).tar.gz
EXTRA_DIST = $(TARFILE)
INCLUDES = $(INCLBENCH)

doit11_SOURCES=doit11.c
doit11_LDADD=sect11.o libdsp.a $(LIBBENCH) @FLIBS@

doit12_SOURCES=doit12.c
doit12_LDADD=sect12.o libdsp.a $(LIBBENCH) @FLIBS@

doit12_fast_SOURCES=doit12-fast.c
doit12_fast_LDADD=sect12.o libdsp.a $(LIBBENCH) @FLIBS@

doit12_ffa_SOURCES=doit12-ffa.c
doit12_ffa_LDADD=sect12.o libdsp.a $(LIBBENCH) @FLIBS@

libdsp_a_SOURCES=sect01.f

libdsp_a_LIBADD=sect04.o

include ../Makefile.common

sect04.f sect11.f sect12.f sect13.f sect14.f sect16.f sect17.f sect18.f sect19.f: extract-mkdir
	@CP@ $(PKG)/$@ $@


