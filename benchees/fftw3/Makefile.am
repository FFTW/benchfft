PRG=doit
PKG=fftw-3.0-attempt2
TARFILE=$(PKG).tar.gz
INCLUDES = $(INCLBENCH) -I$(PKG)/src
LIBFFTW=$(PKG)/src/libfftw.a

if SINGLE
CONFIGURE=./configure --enable-single
else
CONFIGURE=./configure --enable-single
endif

doit_SOURCES=doit.c
doit_LDADD=$(LIBFFTW) $(LIBBENCH) @FLIBS@

doit.o: doit.c extract

$(LIBFFTW): extract
	(cd $(PKG) && $(CONFIGURE) && $(MAKE))

include ../Makefile.common

clean-local:
	rm -rf benchmark extract $(PKG)