BEGIN_BENCH_DOC
BENCH_DOC("name", NAME)
BENCH_DOC("package", "FXT")
BENCH_DOC("version", "26-March-2003")
BENCH_DOC("author", "J&ouml;rg Arndt")
BENCH_DOC("year", "2003")
BENCH_DOC("email", "arndt@jjj.de")
BENCH_DOC("url", "http://www.jjj.de/fxt/")
BENCH_DOC("url-was-valid-on", "Tue Apr 29 18:19:29 EDT 2003")
BENCH_DOC("copyright", 
"Copyright (c) 2001 Joerg Arndt\n"
"\n"
"This program is free software; you can redistribute it and/or modify\n"
"it under the terms of the GNU General Public License as published by\n"
"the Free Software Foundation; either version 2 of the License, or\n"
"(at your option) any later version.\n"
"\n"
"This program is distributed in the hope that it will be useful,\n"
"but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"GNU General Public License for more details.\n"
"\n"
"You should have received a copy of the GNU General Public License\n"
"along with this program; if not, write to the Free Software\n"
"Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA\n")
BENCH_DOC("language", "C++")
BENCH_DOC("notes", "We use whatever CXXFLAGS are selected by fxt's makefile")
BENCH_DOC("notes", NOTES)
#ifdef NO_Complex
BENCH_DOC("notes", "Complex data are stored in separate real/imag arrays.")
#endif
END_BENCH_DOC
