if DOUBLE
PRG=doit
endif

PKG = mixfft05
TARFILE = $(PKG).tar.gz
EXTRA_DIST = $(TARFILE)
LOCAL_COPIES = Mixfft.c
INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
nodist_doit_SOURCES=$(LOCAL_COPIES)
doit_LDADD=$(LIBBENCH) @FLIBS@

$(LOCAL_COPIES): extract
	(echo "#include <config.h>"; tr '\015' ' ' <$(PKG)/$@; echo ""; echo "") | sed -e 's/\/\/.*$$//g' >$@

include ../Makefile.common
