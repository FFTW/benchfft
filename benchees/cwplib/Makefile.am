noinst_PROGRAMS=doit
PKG=cwplib
INCLUDES = -I${top_srcdir}/libbench -I$(PKG)
TARFILE=$(PKG).tar.gz
EXTRA_DIST = $(TARFILE)

doit_SOURCES=doit.c
doit_LDADD=pfafft.o ${top_builddir}/libbench/libbench.a

extract: $(TARFILE)
	gzip -dc ${srcdir}/$(TARFILE) | tar xf -
	@echo timestamp >$@

pfafft.o: $(PKG)/pfafft.c
	$(COMPILE) -c -o $@ $(PKG)/pfafft.c

doit.o: extract
