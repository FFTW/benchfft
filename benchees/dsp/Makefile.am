if SINGLE
if HAVE_F77
PRG=doit11 doit12 doit12-fast doit12-ffa doit14 doit17 doit18
noinst_LIBRARIES=libdsp.a
endif
endif

PKG=dsp
TARFILE=$(PKG).tar.gz
EXTRA_DIST = $(TARFILE) sect17.f sect18.f
INCLUDES = $(INCLBENCH)
LOCAL_COPIES=sect04.f sect11.f sect12.f sect13.f sect14.f sect16.f	\
sect19.f

doit11_SOURCES=doit11.c
doit11_LDADD=sect11.o libdsp.a $(LIBBENCH) @FLIBS@

doit12_SOURCES=doit12.c
doit12_LDADD=sect12.o libdsp.a $(LIBBENCH) @FLIBS@

doit12_fast_SOURCES=doit12-fast.c
doit12_fast_LDADD=sect12.o libdsp.a $(LIBBENCH) @FLIBS@

doit12_ffa_SOURCES=doit12-ffa.c
doit12_ffa_LDADD=sect12.o libdsp.a $(LIBBENCH) @FLIBS@

# Section 1.3: real symmetric transforms.  Not benchmarked.

doit14_SOURCES=doit14.c
doit14_LDADD=sect14.o libdsp.a $(LIBBENCH) @FLIBS@

# Section 1.5: contains PDP11 assembly.  Won't compile unless
# modified. Not benchmarked.

# Section 1.6: Chirp Z transform is outside the scope of the benchmark.

doit17_SOURCES=doit17.c
doit17_LDADD=sect17.o libdsp.a $(LIBBENCH) @FLIBS@

doit18_SOURCES=doit18.c
doit18_LDADD=sect18.o libdsp.a $(LIBBENCH) @FLIBS@


libdsp_a_SOURCES=sect01.f

libdsp_a_LIBADD=sect04.o

include ../Makefile.common

$(LOCAL_COPIES): extract-mkdir
	@CP@ $(PKG)/$@ $@


