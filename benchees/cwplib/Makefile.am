PRG=doit
PKG=cwplib
INCLUDES = $(INCLBENCH) -I$(PKG)
TARFILE=$(PKG).tar.gz
EXTRA_DIST = $(TARFILE)

doit_SOURCES=doit.c
doit_LDADD=pfafft.o $(LIBBENCH)

include ../Makefile.common

pfafft.o: $(PKG)/pfafft.c
	$(COMPILE) -c -o $@ $(PKG)/pfafft.c

$(PKG)/pfafft.c $(PKG)/cwp.h: extract

doit.o: $(PKG)/cwp.h