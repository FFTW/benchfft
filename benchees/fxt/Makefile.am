if DOUBLE
if HAVE_CXX
PRG=doit_dif doit_dit doit_split doit_fht doit_4step doit_twodim doit_ndim
endif
endif

INCLUDES = $(INCLBENCH)
PKG=fxt
TARFILE=$(PKG)-2002.08.26.tar.gz
EXTRA_DIST = $(TARFILE) doit.cc

CPPFLAGS = @CPPFLAGS@ -I$(PKG)/include -I$(PKG)/auxil -I$(PKG)/permute -I$(PKG)/auxbit

LIBFXT = $(PKG)/libfxt.a

doit_dit_SOURCES = doit_dit.cc
doit_dit_LDADD = $(LIBFXT) $(LIBBENCH)
doit_dit.o: doit.cc $(LIBFXT)

doit_dif_SOURCES = doit_dif.cc
doit_dif_LDADD = $(LIBFXT) $(LIBBENCH)
doit_dif.o: doit.cc $(LIBFXT)

doit_split_SOURCES = doit_split.cc
doit_split_LDADD = $(LIBFXT) $(LIBBENCH)
doit_split.o: doit.cc $(LIBFXT)

doit_fht_SOURCES = doit_fht.cc
doit_fht_LDADD = $(LIBFXT) $(LIBBENCH)
doit_fht.o: doit.cc $(LIBFXT)

doit_4step_SOURCES = doit_4step.cc
doit_4step_LDADD = $(LIBFXT) $(LIBBENCH)
doit_4step.o: doit.cc $(LIBFXT)

doit_twodim_SOURCES = doit_twodim.cc
doit_twodim_LDADD = $(LIBFXT) $(LIBBENCH)
doit_twodim.o: doit.cc $(LIBFXT)

doit_ndim_SOURCES = doit_ndim.cc
doit_ndim_LDADD = $(LIBFXT) $(LIBBENCH)
doit_ndim.o: doit.cc $(LIBFXT)

$(LIBFXT): extract
	cd $(PKG) && $(MAKE) CXX="@CXX@" CXXFLAGS="@CXXFLAGS@"

include ../Makefile.common
