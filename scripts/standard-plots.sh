#! /bin/sh

data=$1
dname=`basename $data .speed`

for rank in 1 2 3; do
    for problem in dcxx drxx scxx scxx srxx; do
	for p2 in p2 np2; do
	    echo "outputting ${dname}.${rank}d.${problem}.${p2}.ps"
	    pat=`echo $problem | sed 's/xy/[io][fb]/g'`
	    (sh plot-title.sh $problem $rank $p2;
		egrep "$pat" $data | perl grep-rank.pl $rank |
		perl grep-${p2}.pl | perl grace-plot.pl --no-dups) |
		gracebat -pipe -printfile ${dname}.${rank}d.${problem}.${p2}.ps
	done
    done
done
