if DOUBLE
if HAVE_CXX
PRG = doit
endif
endif

INCLUDES = $(INCLBENCH) -I$(PKG)/include
PKG = arprec
TARFILE = $(PKG).tar.gz
EXTRA_DIST = $(TARFILE) mp.h.patch mpreal.h.patch

ARPREC = $(PKG)/lib/libmp.a

doit_SOURCES=doit.cc $(FFT)
doit_LDADD=$(ARPREC) $(LIBBENCH)
doit.o: $(ARPREC)

$(ARPREC): extract-and-patch
	cd $(PKG) && ./configure && $(MAKE)

extract-and-patch: $(TARFILE)
	(@GZIP@ -dc ${srcdir}/$(TARFILE) | @TAR@ xf -) && echo timestamp >$@
	$(PATCH) $(PKG)/include/mp/mp.h < $(srcdir)/mp.h.patch
	$(PATCH) $(PKG)/include/mp/mpreal.h < $(srcdir)/mpreal.h.patch

include ../Makefile.common

clean-local::
	rm -f extract-and-patch