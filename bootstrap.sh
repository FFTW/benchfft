# script to initialize automake/autoconf etc
echo "Please ignore warnings and errors"
autoheader2.50
aclocal
automake --add-missing
autoconf2.50
autoheader2.50
aclocal
automake --add-missing
autoconf2.50

rm -f config.cache

# Configure in separate directory so as not to mess source code
(
    rm -rf OBJ
    mkdir OBJ
    cd OBJ
    ../configure $*
)
