#! /bin/sh

# Usage: desc2web.sh <foo.desc> [<directory>]
# Creates a directory foo in <directory> (or .) containing an index.html
# and images etcetera for the data foo.tar.gz, described by foo.desc.

sd=`dirname $0`
case $sd in
    /*) ;;
    *) sd=`pwd`/$sd ;;
esac

desc=$1
tgz=`echo $desc | sed 's/\.desc/.tar.gz/'`
dir=${2-`pwd`}
name=`basename $desc .desc`
type=${3-speed}

abbrevname=`head -1 $desc | cut -d: -f1`

rm -rf $dir/$name $dir/$abbrevname
mkdir $dir/$name || exit 1
mkdir $dir/$name/data
ln -s $dir/$name $dir/$abbrevname

cp $tgz $dir/$name/data

cat > $dir/$name/index.html <<EOF
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
     "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
EOF

echo '<title>' >> $dir/$name/index.html
head -1 $desc |sed 's/^.*: *//' >> $dir/$name/index.html
echo '</title>' >> $dir/$name/index.html

# fixme: add bookmark links, favicon, "go back" link, etcetera
cat >> $dir/$name/index.html <<EOF
</head>
<body text="#000000" bgcolor="#ffffff">
EOF

echo '<h1 align="center">' >> $dir/$name/index.html
head -1 $desc |sed 's/^.*: *//' >> $dir/$name/index.html
echo '</h1>' >> $dir/$name/index.html

echo "<p>" >> $dir/$name/index.html
perl -p -e 's/\n/<p>/g' $desc | perl -p -e 's/<p>(<p>)?(<p>)*/\n\1/g' |
    tail +2 | perl -p -e 's/<p>/\n<p>/g' >> $dir/$name/index.html

cd $dir/$name/data
tgz=`basename $tgz`
tar xzf $tgz

# sanitize
verboten='djbfft|athfft|pfftw'
for suff in speed accuracy config.h config.log config.status; do
    egrep -i -v "$verboten" ${name}.${suff} > foo && mv -f foo ${name}.${suff}
done
guile -s $sd/sanitize.scm < ${name}.info > foo && mv -f foo ${name}.info
rm $tgz
tar czf $tgz ${name}.*

cd ..

echo '<p>Compilers and flags (unless overridden):</p> <ul>' >> index.html

echo -n '<li>C: <code>' >> index.html
grep "(cc" data/${name}.info |head -1 |sed 's/^(cc *\"//;s/\")$//' |tr -d '\n' >> index.html
echo '</code>' >> index.html

echo -n '<li>C++: <code>' >> index.html
grep "(cxx" data/${name}.info |head -1 |sed 's/^(cxx *\"//;s/\")$//' |tr -d '\n' >> index.html
echo '</code>' >> index.html

echo -n '<li>Fortran: <code>' >> index.html
grep "(f77" data/${name}.info |head -1 |sed 's/^(f77 *\"//;s/\")$//' |tr -d '\n' >> index.html
echo '</code>' >> index.html

echo '</ul>' >> index.html

echo '<p>Raw data files: <a href="ftp://ftp.fftw.org/pub/fftw/numbers/'$tgz'">'$tgz'</a>' >> index.html

echo '<hr>' >> index.html

sh $sd/standard-plots.sh data/${name}.${type} |while read ps; do
    png=`basename $ps .ps`.png
    convert -rotate 90 $ps $png
    rm $ps
    width=`identify $png |cut -d' ' -f3 |cut -dx -f1`
    height=`identify $png |cut -d' ' -f3 |cut -dx -f2`
    echo '<p align="center"><img src="'$png'" width="'$width'" height="'$height'">' >> index.html
done

cp -f data/$tgz /home/fftw/ftp/numbers/$tgz
rm -rf data

cat >> index.html <<EOF
</body>
</html>
EOF
