if SINGLE
if HAVE_F77
PRG=doit
noinst_LIBRARIES=libfft.a
endif
endif

EXTRA_DIST = morris.f morris.f.orig
INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=libfft.a $(LIBBENCH) @FLIBS@

libfft_a_SOURCES = morris.f

include ../Makefile.common

$(LOCAL_COPIES): extract-mkdir
	@CP@ $(PKG)/$@ $@


