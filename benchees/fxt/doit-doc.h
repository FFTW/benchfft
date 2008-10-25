BEGIN_BENCH_DOC
BENCH_DOC("name", NAME)
BENCH_DOC("package", "FXT")
BENCH_DOC("version", "31-August-2008")
BENCH_DOC("author", "J&ouml;rg Arndt")
BENCH_DOC("year", "2008")
BENCH_DOC("email", "arndt@jjj.de")
BENCH_DOC("url", "http://www.jjj.de/fxt/")
BENCH_DOC("url-was-valid-on", "Sat Oct 25 15:56:44 EDT 2008")
BENCH_DOC("copyright", 
	  "fxt is distributed under the GNU GENERAL PUBLIC LICENSE (GPL)\n"
	  "cf. the file COPYING.txt")
BENCH_DOC("language", "C++")
BENCH_DOC("notes", "For GNU g++, we use whatever CXXFLAGS are selected by fxt's makefile")
BENCH_DOC("notes", NOTES)
#ifdef NO_Complex
BENCH_DOC("notes", "Complex data are stored in separate real/imag arrays.")
#endif
END_BENCH_DOC
