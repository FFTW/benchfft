#! /bin/sh

sd=`dirname $0`
case $sd in
    /*) ;;
    *) sd=`pwd`/$sd ;;
esac

data=$1
case $data in
    *.speed) type=speed; dname=`basename $data .speed` ;;
    *.accuracy) type=accuracy; dname=`basename $data .accuracy` ;;
esac

if test $type = speed; then

plotflags="--no-dups"

for rank in 1 2 3; do
    for problem in dcxx drxx; do
	for p2 in p2 np2; do
	    pat=`echo $problem | sed 's/xx/[io][fb]/g'`
	    (sh $sd/plot-title.sh $problem $rank $p2;
		egrep "$pat" $data |
		perl $sd/grep-rank.pl $rank |
		perl $sd/grep-${p2}.pl | perl $sd/grace-plot.pl $plotflags) |
		gracebat -pipe -printfile ${dname}.${rank}d.${problem}.${p2}.ps
	    echo "${dname}.${rank}d.${problem}.${p2}.ps"
	done
    done
done

for rank in 1 2 3; do
    for problem in scxx srxx; do
	for p2 in p2 np2; do
	    pat=`echo $problem | sed 's/xx/[io][fb]/g'`
	    (sh $sd/plot-title.sh $problem $rank $p2;
		egrep "$pat" $data |
		perl $sd/grep-rank.pl $rank |
		perl $sd/grep-${p2}.pl | perl $sd/grace-plot.pl $plotflags) |
		gracebat -pipe -printfile ${dname}.${rank}d.${problem}.${p2}.ps
	    echo "${dname}.${rank}d.${problem}.${p2}.ps"
	done
    done
done

else # $type = accuracy

rank=1

plotflags="--no-dups --plot-worst --accuracy"

for problem in dcxx drxx scxx srxx; do
    pat=`echo $problem | sed 's/xx/[io][fb]/g'`

    test `egrep "$pat" $data | perl $sd/grep-rank.pl $rank | wc -l` -gt 1 || continue

    (sh $sd/plot-title.sh $problem $rank;
    egrep "$pat" $data |
	perl $sd/grep-rank.pl $rank |
	perl $sd/grace-plot.pl $plotflags) |
	gracebat -pipe -printfile ${dname}.${rank}d.${problem}.ps
    echo "${dname}.${rank}d.${problem}.ps"

    (sh $sd/plot-title.sh $problem $rank;
     echo "@subtitle \"enlargement of most accurate FFTs\"";
     egrep "$pat" $data |
	perl $sd/grep-rank.pl $rank |
	perl $sd/grace-plot.pl $plotflags --accurate-only) |
	gracebat -pipe -printfile ${dname}.${rank}d.${problem}.acc.ps
    echo "${dname}.${rank}d.${problem}.acc.ps"

    (sh $sd/plot-title.sh $problem $rank;
     echo "@subtitle \"enlargement of most accurate FFTs, powers of two only\"";
     egrep "$pat" $data |
	perl $sd/grep-rank.pl $rank | perl $sd/grep-p2.pl |
	perl $sd/grace-plot.pl $plotflags --accurate-only) |
	gracebat -pipe -printfile ${dname}.${rank}d.${problem}.acc.p2.ps
    echo "${dname}.${rank}d.${problem}.acc.p2.ps"
done

fi
