if DOUBLE
PRG=doit0 doit1 doit2 doit3 doit4 doit5 doit7
endif

INCLUDES = $(INCLBENCH)
PKG=fft
TARFILE=$(PKG)_source.tar.gz
EXTRA_DIST = $(TARFILE) doit.c

doit0_SOURCES = doit0.c
doit0_LDADD = $(PKG)/mayer/fft.o $(PKG)/mayer/fht0.o $(LIBBENCH)
doit0.o: doit.c

doit1_SOURCES = doit1.c
doit1_LDADD = $(PKG)/mayer/fft.o $(PKG)/mayer/fht1.o $(LIBBENCH)
doit1.o: doit.c

doit2_SOURCES = doit2.c
doit2_LDADD = $(PKG)/mayer/fft.o $(PKG)/mayer/fht2.o $(LIBBENCH)
doit2.o: doit.c

doit3_SOURCES = doit3.c
doit3_LDADD = $(PKG)/mayer/fft.o $(PKG)/mayer/fht3.o $(LIBBENCH)
doit3.o: doit.c

doit4_SOURCES = doit4.c
doit4_LDADD = $(PKG)/mayer/fft.o $(PKG)/mayer/fht4.o $(LIBBENCH)
doit4.o: doit.c

doit5_SOURCES = doit5.c
doit5_LDADD = $(PKG)/mayer/fft.o $(PKG)/mayer/fht5.o $(LIBBENCH)
doit5.o: doit.c

doit7_SOURCES = doit7.c
doit7_LDADD = $(PKG)/mayer/fft.o $(PKG)/mayer/fht7.o $(LIBBENCH)
doit7.o: doit.c

MAYER_OBJS = $(PKG)/mayer/fft.o $(PKG)/mayer/fht0.o		\
$(PKG)/mayer/fht1.o $(PKG)/mayer/fht2.o $(PKG)/mayer/fht2d.o	\
$(PKG)/mayer/fht3.o $(PKG)/mayer/fht4.o $(PKG)/mayer/fht5.o	\
$(PKG)/mayer/fht6.o $(PKG)/mayer/fht7.o

$(MAYER_OBJS): extract
	cd $(PKG)/mayer && $(MAKE) CC="@CC@" CFLAGS="@CFLAGS@"

include ../Makefile.common
