if DOUBLE
if HAVE_NAG
PRG = doit_ws doit_ip doit_mult doit_cmplx
endif
endif

EXTRA_DIST = doit.c

INCLUDES = $(INCLBENCH)

doit_ws_SOURCES = doit-ws.c
doit_ws_LDADD = @LIBNAG@ $(LIBBENCH) @FLIBS@
doit-ws.o: doit.c

doit_ip_SOURCES = doit-ip.c
doit_ip_LDADD = @LIBNAG@ $(LIBBENCH) @FLIBS@
doit-ip.o: doit.c

doit_mult_SOURCES = doit-mult.c
doit_mult_LDADD = @LIBNAG@ $(LIBBENCH) @FLIBS@
doit-mult.o: doit.c

doit_cmplx_SOURCES = doit-cmplx.c
doit_cmplx_LDADD = @LIBNAG@ $(LIBBENCH) @FLIBS@

include ../Makefile.common
