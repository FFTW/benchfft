if SINGLE
PRG=doit
endif

PKG=vBigDSP
TARFILE=$(PKG)-20030322.tgz
EXTRA_DIST = $(TARFILE)

INCLUDES = $(INCLBENCH) -I$(PKG)

doit_SOURCES=doit.c
nodist_doit_SOURCES = $(PKG)/vBigDSP.c $(PKG)/vBigDSP.h
doit_CFLAGS=$(AM_CFLAGS) @ALTIVEC_CFLAGS@ -Dvalloc=bench_malloc
doit_LDADD=$(LIBBENCH)

# I don't understand this: automake produces the name doit-doit.o
# internally
doit-doit.o : extract

include ../Makefile.common
