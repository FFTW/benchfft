#! /bin/sh

echo -n "@title \""

if test $# -gt 3; then
    echo -n "$4"
fi

problem=$1
case $problem in
    dci*) echo -n "double-precision complex, in-place" ;;
    dco*) echo -n "double-precision complex, out-of-place" ;;
    sci*) echo -n "single-precision complex, in-place" ;;
    sco*) echo -n "single-precision complex, out-of-place" ;;
    dri*) echo -n "double-precision real-data, in-place" ;;
    dro*) echo -n "double-precision real-data, out-of-place" ;;
    sri*) echo -n "single-precision real-data, in-place" ;;
    sro*) echo -n "single-precision real-data, out-of-place" ;;
    dcx*) echo -n "double-precision complex" ;;
    scx*) echo -n "single-precision complex" ;;
    drx*) echo -n "double-precision real-data" ;;
    srx*) echo -n "single-precision real-data" ;;
esac

rank=$2
echo -n ", ${rank}d transforms"

echo "\""

if test $# -gt 2; then
    p2=$3
    case $p2 in
	p2) echo "@subtitle \"powers of two\"" ;;
	np2) echo "@subtitle \"non-powers of two\"" ;;
    esac
fi
