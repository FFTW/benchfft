#! /bin/sh

data=$1
dname=`basename $data .speed`

for problem in dcxf dcif dcof drxf drif drof scxf scif scof srxf srif srof; do
    for p2 in p2 np2; do
	echo "outputting ${dname}.1d.${problem}.${p2}.ps"
	pat=`echo $problem | sed 's/x/[io]/g'`
	(sh plot-title.sh $problem 1 $p2;
	    egrep "$pat" $data | perl grep-rank.pl 1 |
	    perl grep-${p2}.pl | perl grace-plot.pl --no-dups) |
	    gracebat -pipe -printfile ${dname}.1d.${problem}.${p2}.ps
    done
done

for rank in 2 3; do
    for problem in dcif drif scif srif; do
	for p2 in p2 np2; do
	    echo "outputting ${dname}.${rank}d.${problem}.${p2}.ps"
	    (sh plot-title.sh $problem $rank $p2;
		egrep $problem $data | perl grep-rank.pl $rank |
		perl grep-${p2}.pl | perl grace-plot.pl --no-dups) |
		gracebat -pipe -printfile ${dname}.${rank}d.${problem}.${p2}.ps
	done
    done
done
