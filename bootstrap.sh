# script to initialize automake/autoconf etc
echo "PLEASE IGNORE WARNINGS AND ERRORS"

# paranoia: sometimes autoconf doesn't get things right the first time
rm -rf autom4te.cache
autoreconf --verbose --install --symlink --force

rm -f config.cache

# Configure in separate directory so as not to mess source code
(
    rm -rf OBJ
    mkdir OBJ
    cd OBJ
    ../configure $*
)
