if DOUBLE
if HAVE_F77
PRG=doit
endif
endif

PKG = f77-ieee
TARFILE = $(PKG).tar.gz
EXTRA_DIST = $(TARFILE) mpfun.patch
LOCAL_COPIES = mpfun.f

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
nodist_doit_SOURCES=mpfun.f
doit_LDADD=$(LIBBENCH) @FLIBS@

mpfun.f: extract
	@CP@ -f $(PKG)/$@ $@
	@PATCH@ $@ $(srcdir)/mpfun.patch

include ../Makefile.common
