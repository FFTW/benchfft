if HAVE_F77
PRG=doit
endif

PKG = f77-ieee
TARFILE = $(PKG).tar.gz
EXTRA_DIST = $(TARFILE)
LOCAL_COPIES = mpfun.f

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c mpfun.f
doit_LDADD=$(LIBBENCH) @FLIBS@

mpfun.f: extract $(CONFIG_HEADER)
	cp -f $(PKG)/mpfun.f $@
	patch $@ $(srcdir)/mpfun.patch

include ../Makefile.common
