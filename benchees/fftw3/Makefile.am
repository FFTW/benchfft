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

CC_FOR_FFTW3=@CC@

doit_SOURCES=doit.c
doit_LDADD=$(LIBFFTW) $(LIBBENCH) @FLIBS@

doit.o: doit.c $(LIBFFTW)

$(LIBFFTW): extract
	(cd $(PKG) && (CC="$(CC_FOR_FFTW3)" $(CONFIGURE)) && $(MAKE))

include ../Makefile.common

clean-local:
	rm -rf benchmark extract $(PKG) *.o $(PRG)