if HAVE_FFMPEG
if SINGLE
PRG=doit
endif
endif

INCLUDES = $(INCLBENCH) `@FFMPEG_CONFIG@ --cflags`

doit_SOURCES=doit.c
doit_LDADD=$(LIBBENCH) $(FFMPEG_LIBS)

include ../Makefile.common
