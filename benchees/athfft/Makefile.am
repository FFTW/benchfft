PRG=doit
PKG=athfft-0.99
TARFILE=$(PKG).tar.gz
INCLUDES = $(INCLBENCH) -I$(PKG)
LIBFFT=$(PKG)/libathfft.a
EXTRA_DIST = $(TARFILE)

if SINGLE
CONFIGURE=./configure --enable-single
else
CONFIGURE=./configure
endif

CC_FOR_ATHFFT=@CC@

doit_SOURCES=doit.c
doit_LDADD=$(LIBFFT) $(LIBBENCH)

doit.o: doit.c $(LIBFFT)

$(LIBFFT): extract
	(cd $(PKG) && (CC="$(CC_FOR_ATHFFT)" $(CONFIGURE)) && $(MAKE))

include ../Makefile.common

clean-local::
	rm -rf benchmark extract $(PKG) *.o $(PRG)
