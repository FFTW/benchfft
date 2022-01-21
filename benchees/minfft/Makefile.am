PRG=doit

AM_CPPFLAGS=$(INCLBENCH)

if SINGLE
AM_CFLAGS=-DMINFFT_SINGLE
endif

doit_SOURCES=doit.c minfft.c minfft.h
doit_LDADD=$(LIBBENCH)

include ../Makefile.common
