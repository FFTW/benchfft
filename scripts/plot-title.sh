#! /bin/sh

printf "@title \""

if test $# -gt 3; then
    printf "$4"
fi

problem=$1
case $problem in
    dci*) printf "double-precision complex, in-place" ;;
    dco*) printf "double-precision complex, out-of-place" ;;
    sci*) printf "single-precision complex, in-place" ;;
    sco*) printf "single-precision complex, out-of-place" ;;
    dri*) printf "double-precision real-data, in-place" ;;
    dro*) printf "double-precision real-data, out-of-place" ;;
    sri*) printf "single-precision real-data, in-place" ;;
    sro*) printf "single-precision real-data, out-of-place" ;;
    dcx*) printf "double-precision complex" ;;
    scx*) printf "single-precision complex" ;;
    drx*) printf "double-precision real-data" ;;
    srx*) printf "single-precision real-data" ;;
esac

rank=$2
printf ", ${rank}d transforms"

echo "\""

if test $# -gt 2; then
    p2=$3
    case $p2 in
	p2) echo "@subtitle \"powers of two\"" ;;
	np2) echo "@subtitle \"non-powers of two\"" ;;
    esac
fi
