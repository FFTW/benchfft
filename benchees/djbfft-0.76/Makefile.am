noinst_PROGRAMS=doit
PKG=djbfft-0.76
INCLUDES = -I${top_srcdir}/libbench -I$(PKG)
TARFILE=$(PKG).tar.gz
EXTRA_DIST = $(TARFILE)

doit_SOURCES=doit.c
doit_LDADD=$(PKG)/djbfft.a ${top_builddir}/libbench/libbench.a

doit.o: doit.c extract $(CONFIG_HEADER) $(PKG)/djbfft.a
	$(COMPILE) \
	-DCC="\"`cat $(PKG)/conf-cc`\"" \
	-DCONF_OPT="\"`cat $(PKG)/conf-opt`\"" \
	-DAUTO_OPT="\"`cat $(PKG)/auto_opt`\"" \
         -c $<

$(PKG)/djbfft.a: $(PKG)
	(cd $(PKG); make)

extract: $(TARFILE)
	gzip -dc ${srcdir}/$(TARFILE) | tar xf -
	@echo timestamp >$@

$(PKG)/fftc4.h $(PKG)/fftc8.h $(PKG)/fftfreq.h: $(PKG)