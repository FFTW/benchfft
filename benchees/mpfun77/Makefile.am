if HAVE_F77
PRG=doit
endif

PKG = f77-ieee
TARFILE = $(PKG).tar.gz
EXTRA_DIST = $(TARFILE) mpfun.patch
LOCAL_COPIES = mpfun.f

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c mpfun.f
doit_LDADD=$(LIBBENCH) @FLIBS@

mpfun.f: extract
	cp -f $(PKG)/$@ $@
	patch $@ $(srcdir)/mpfun.patch

include ../Makefile.common
