PRG=doit
PKG=fftw-3.0-attempt2
TARFILE=$(PKG).tar.gz
INCLUDES = $(INCLBENCH) -I$(PKG)/src
LIBFFTW=$(PKG)/src/libfftw.a

if SINGLE
CONFIGURE=./configure --enable-single
else
CONFIGURE=./configure
endif

doit_SOURCES=doit.c
doit_LDADD=$(LIBFFTW) $(LIBBENCH) @FLIBS@

doit.o: doit.c extract

$(LIBFFTW): extract
	(cd $(PKG) && (CC=$(CC) CFLAGS=$(CFLAGS) $(CONFIGURE)) && $(MAKE))

include ../Makefile.common

clean-local:
	rm -rf benchmark extract $(PKG) *.o $(PRG)