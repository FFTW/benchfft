if HAVE_PFFTW_ARCH
PRG=doit
endif

PKG=pfftw-0.03
TARFILE=$(PKG).tar.gz
PKGARCH=$(PKG)/@PFFTW_ARCH@
AM_CPPFLAGS = $(INCLBENCH) -I$(PKGARCH)
LIBPFFTW=libpfftw.a
LIBFFT=$(PKGARCH)/$(LIBPFFTW)
EXTRA_DIST = $(TARFILE)

doit_SOURCES=doit.c
doit_LDADD=$(LIBFFT) $(LIBBENCH)

doit.o: doit.c $(LIBFFT)

$(LIBFFT): extract
	(cd $(PKGARCH) && $(MAKE) $(LIBPFFTW))

include ../Makefile.common

clean-local::
	rm -rf benchmark extract $(PKG) *.o $(PRG)
