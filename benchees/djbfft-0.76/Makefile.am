PRG=doit
PKG=djbfft-0.76
INCLUDES = $(INCLBENCH) -I$(PKG)
TARFILE=$(PKG).tar.gz
# EXTRA_DIST = $(TARFILE)

doit_SOURCES=doit.c
doit_LDADD=$(PKG)/djbfft.a $(LIBBENCH)

include ../Makefile.common

doit.o: doit.c extract $(CONFIG_HEADER) $(PKG)/djbfft.a
	$(COMPILE) \
	-DDJBFFT_CC="\"`cat $(PKG)/conf-cc`\"" \
	-DCONF_OPT="\"`cat $(PKG)/conf-opt`\"" \
	-DAUTO_OPT="\"`cat $(PKG)/auto_opt`\"" \
         -c $<

$(PKG)/djbfft.a: extract
	(cd $(PKG); $(MAKE) djbfft.a)

$(PKG)/fftc4.h $(PKG)/fftc8.h $(PKG)/fftfreq.h: extract