if HAVE_F77
PRG=doit11 doit12
noinst_LIBRARIES=libdsp.a
endif

PKG=dsp
TARFILE=$(PKG).tar.gz
EXTRA_DIST = $(TARFILE)
INCLUDES = $(INCLBENCH)

doit11_SOURCES=doit11.c
doit11_LDADD=libdsp.a $(LIBBENCH) @FLIBS@

doit12_SOURCES=doit12.c
doit12_LDADD=libdsp.a $(LIBBENCH) @FLIBS@

libdsp_a_SOURCES=sect01.f

# sect15.f is Fraser's program which does not compile as distributed
# (It even contains PDP assembly)
libdsp_a_LIBADD=sect04.o sect11.o sect12.o sect13.o sect14.o sect16.o	\
sect17.o sect18.o sect19.o

include ../Makefile.common

sect04.f sect11.f sect12.f sect13.f sect14.f sect16.f sect17.f sect18.f sect19.f: extract-mkdir
	@CP@ $(PKG)/$@ $@


