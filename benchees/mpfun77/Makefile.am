if DOUBLE
if HAVE_F77
PRG=doit
noinst_LIBRARIES=libfft.a
endif
endif

PKG = f77-ieee
TARFILE = $(PKG).tar.gz
EXTRA_DIST = $(TARFILE) mpfun.patch
LOCAL_COPIES = mpfun.f

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=libfft.a $(LIBBENCH) @FLIBS@

nodist_libfft_a_SOURCES=mpfun.f

mpfun.f: extract
	@CP@ -f $(PKG)/$@ $@
	@PATCH@ $@ <$(srcdir)/mpfun.patch

include ../Makefile.common
