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

rm -rf $dir/$name
mkdir $dir/$name || exit 1

cp $tgz $dir/$name

cat > $dir/$name/index.html <<EOF
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
     "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
EOF

echo '<title>' >> $dir/$name/index.html
head -1 $desc >> $dir/$name/index.html
echo '</title>' >> $dir/$name/index.html

# fixme: add bookmark links, favicon, "go back" link, etcetera
cat >> $dir/$name/index.html <<EOF
</head>
<body text="#000000" bgcolor="#ffffff">
EOF

echo '<h1 align="center">' >> $dir/$name/index.html
head -1 $desc >> $dir/$name/index.html
echo '</h1>' >> $dir/$name/index.html

echo "<p>" >> $dir/$name/index.html
tail +2 $desc >> $dir/$name/index.html

cd $dir/$name
tgz=`basename $tgz`

tar xzf $tgz
mkdir data
mv `ls $name.* |grep -v $tgz` data

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

cp -f $tgz /home/fftw/ftp/numbers/$tgz
rm -rf data $tgz

cat >> index.html <<EOF
</body>
</html>
EOF
