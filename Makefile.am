AUTOMAKE_OPTIONS = gnu

SUBDIRS = libbench benchees scripts

EXTRA_DIST = COPYRIGHT bootstrap.sh

collect:
	$(SHELL) ${top_srcdir}/scripts/collect ${top_builddir}/benchees

# this is a hack
RECURSIVE_TARGETS =
RECURSIVE_TARGETS += benchmark-recursive accuracy-recursive
benchmark: benchmark-recursive
benchmark-am:
accuracy: accuracy-recursive
accuracy-am:

