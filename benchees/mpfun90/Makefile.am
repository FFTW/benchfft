if DOUBLE
if HAVE_F90
PRG = doit
endif
endif

PKG = f90-ieee
TARFILE = $(PKG).tar.gz
EXTRA_DIST = $(TARFILE) mpfun90.patch
LOCAL_COPIES = mpfun90.f90

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=mpfun90.o $(LIBBENCH) @F90LIBS@

# We need to extract just modules mpfuna and mpfung, and "unmodulize"
# mpfung.  For one thing, autoconf doesn't know how to mangle
# identifiers for F90 modules.  Also VAST/f90 has problems compiling
# some of the extraneous routines.  (Sigh.)

mpfun90.f90: extract
	rm -f $@
	head -133 $(PKG)/mpfun90.f > $@
	tail +6081 $(PKG)/mpfun90.f | head -630 >> $@
	@PATCH@ $@ <$(srcdir)/mpfun90.patch

mpfun90.o: mpfun90.f90
	$(F90) @F90FLAGS@ -c mpfun90.f90

include ../Makefile.common
