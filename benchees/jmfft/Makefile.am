if DOUBLE
PRG=doit
endif

PKG = JMFFTC-1.0
TARFILE = jmfftc-1.0.tar.gz
EXTRA_DIST = $(TARFILE) jmfft-7.2.tar.gz

INCLUDES = $(INCLBENCH) -I$(PKG)/Source

LIBJMFFT = $(PKG)/Source/libjmfft.a

doit_SOURCES = doit.c
doit_LDADD = $(LIBJMFFT) $(LIBBENCH) @FLIBS@
doit.o: extract

$(LIBJMFFT): extract
	cd $(PKG)/Source && $(MAKE) CC="$(CC)" CFLAGS="-c $(CFLAGS)" AR="$(AR)" ARFLAGS="cru" RM="$(RM)" RMFLAGS="-rf" && $(RANLIB) libjmfft.a

include ../Makefile.common
