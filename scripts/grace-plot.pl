#! /usr/bin/perl

$no_dups = 0;
$accuracy = 0;
$accurate_only = 0;
while (@ARGV) {
    $arg = shift;
    $no_dups = 1 if ($arg eq "--no-dups");
    $accuracy = 1 if ($arg eq "--accuracy");
    $accurate_only = 1 if ($arg eq "--accurate-only");
}

#############################################################################

# Here, we make a map from transform names to line styles, so that each
# transform is plotted in a consistent style.  Which style we assign to
# which FFT is fairly arbitrary, but we should try to make FFTs that are
# close in performance visually distinct on the graph.  We also make
# FFTs of related heritage the same color.  Jokes wherever possible:
# bloodworth is red, green is green, cross is a plus, etcetera.

# A style consists of (: separated):
#   line-color:line-style:line-width
#   :symbol-color:symbol-shape:symbol-size:symbol-fill-color

# The colors can be any one of:
#    none, white, black, red, green, blue, yellow, brown, grey, violet,
#          cyan, magenta, orange, indigo, maroon, turquoise, green4

# The line-style can be any one of:
#    none, solid, dot, dash, dash2, dotdash, dotdash2, dotdotdash, dotdashdash

# The symbol-shape can be any one of:
#    none, circle, square, diamond, triangle-up, triangle-down, 
#    triangle-right, plus, x, star, '<char>'
# ('<char>' indicates that the letter given by <char> is used as the symbol)

# (The default symbol-size/line-width is 1.0.)

%styles = (
	   "arprec" => "green4:solid:1:green4:circle:0.5:none",
	   "athfft" => "indigo:dot:1:indigo:circle:0.5:none",
	   "bloodworth" => "red:dot:1:red:square:0.5:none",
	   "bloodworth-fht" => "red:solid:1:red:square:0.5:red",
	   "burrus-sffteu" => "green:solid:1:green:circle:0.5:none",
	   "cross" => "black:solid:1:black:plus:0.5:none",
	   "cwplib" => "orange:dotdash:1:orange:square:0.5:none",
	   "dfftpack" => "brown:dash:1:brown:triangle-right:0.5:none",
	   "djbfft-0.76" => "maroon:solid:1:maroon:square:0.5:none",
	   "dsp79-rader" => "yellow:solid:1:yellow:circle:0.5:none",
	   "dsp79-FAST" => "yellow:dot:1:yellow:circle:0.5:yellow",
	   "dsp79-FFA" => "yellow:solid:1:yellow:x:0.5:none",
	   "dsp79-FFT842" => "yellow:dot:1:yellow:triangle-up:0.5:none",
	   "dsp79-singleton" => "yellow:solid:1:yellow:star:0.8:none",
	   "dsp79-wfta" => "yellow:dash:1:yellow:square:0.5:none",
	   "dsp79-morris" => "yellow:dotdash:1:yellow:square:0.5:yellow",
	   "dxml" => "black:solid:1:black:star:0.5:none",
	   "emayer" => "turquoise:dot:1.5:turquoise:plus:0.5:none",
	   "esrfft" => "indigo:dash:1:indigo:square:0.4:indigo",
	   "essl" => "black:solid:1:black:star:0.5:none",
	   "ffte" => "indigo:solid:2:indigo:none:0.5:none",
	   "fftpack" => "brown:dash:1:brown:triangle-right:0.5:brown",
	   "fftreal" => "red:dash:1:red:circle:0.5:none",
	   "fftw2-measure" => "cyan:solid:1:cyan:circle:0.5:cyan",
	   "fftw2-estimate" => "cyan:solid:1:cyan:circle:0.5:none",
	   "fftw2-nd-measure" => "cyan:dash:1:cyan:circle:0.5:cyan",
	   "fftw2-nd-estimate" => "cyan:dash:1:cyan:circle:0.5:none",
	   "fftw3_simd" => "black:solid:1:blue:circle:0.5:blue",
	   "fftw3_simd-impatient" => "black:solid:1:blue:triangle-up:0.5:none",
	   "fftw3_simd-estimate" => "black:solid:1:blue:circle:0.5:none",
	   "fftw3" => "blue:solid:1:blue:circle:0.5:blue",
	   "fftw3-impatient" => "blue:solid:1:blue:triangle-up:0.5:none",
	   "fftw3-estimate" => "blue:solid:1:blue:circle:0.5:none",
	   "fftw3 out-of-place" => "blue:solid:1:blue:circle:0.5:blue",
	   "fftw3-impatient out-of-place" => "blue:solid:1:blue:triangle-up:0.5:none",
	   "fftw3-estimate out-of-place" => "blue:solid:1:blue:circle:0.5:none",
	   "fftw3 in-place" => "blue:dot:1:blue:circle:0.3:blue",
	   "fftw3-impatient in-place" => "blue:dot:1:blue:triangle-up:0.3:none",
	   "fftw3-estimate in-place" => "blue:dot:1:blue:circle:0.3:none",
	   "fftw3-r2r in-place" => "blue:dotdash:1:blue:square:0.4:none",
	   "fftw3-r2r out-of-place" => "blue:dash:1:blue:square:0.5:blue",
	   "fftw3-r2r" => "blue:dash:1:blue:square:0.5:blue",
	   "fxt-4step" => "yellow:solid:2:yellow:circle:0.25:none",
	   "fxt-matrixfft" => "yellow:solid:2:yellow:circle:0.25:none",
	   "fxt-dif" => "yellow:solid:1:yellow:circle:0.5:none",
	   "fxt-dit" => "yellow:solid:1:yellow:circle:0.5:yellow",
	   "fxt-fht" => "yellow:dot:1:yellow:triangle-up:0.5:none",
	   "fxt-ndim" => "yellow:dash:1:yellow:star:1:none",
	   "fxt-split" => "yellow:dotdotdash:1:yellow:square:0.7:none",
	   "fxt-twodim" => "yellow:dot:1:yellow:diamond:0.7:none",
	   "goedecker" => "brown:solid:2:brown:circle:0.25:none",
	   "gpfa" => "maroon:dash:1:maroon:diamond:0.5:none",
	   "gpfa-3d" => "maroon:dash:1:maroon:diamond:0.5:none",
	   "green-ffts-2.0" => "green:solid:2:green:circle:0.25:green",
	   "gsl-mixed-radix" => "violet:solid:2:violet:circle:0.25:none",
	   "gsl-radix2" => "violet:dash:1:violet:circle:0.5:violet",
	   "gsl-radix2-dif" => "violet:dash:1:violet:star:0.5:none",
	   "gsl-radix2-dit" => "violet:dot:1:violet:circle:0.5:none",
	   "harm" => "green4:dot:1:green4:star:0.5:none",
	   "intel-ipps" => "black:solid:2:black:circle:0.2:black",
	   "intel-mkl" => "black:solid:1:black:circle:0.5:black",
	   "intel-mkl-dfti" => "black:solid:1:black:square:0.5:none",
	   "intel-mkl-f" => "black:dash:1:black:circle:0.5:none",
	   "jmfftc" => "turquoise:solid:2:turquoise:circle:0.5:turquoise",
	   "krukar" => "cyan:dash:1:cyan:triangle-up:0.5:none",
	   "mfft" => "orange:solid:2:orange:circle:0.25:none",
	   "monnier" => "turquoise:dash:1:turquoise:square:0.5:none",
	   "morris82" => "turquoise:solid:1:turquoise:star:0.4:blue",
	   "mixfft" => "blue:dot:1:blue:star:0.5:none",
	   "bailey" => "green4:dot:1:green4:square:0.5:green4",
	   "mpfun77" => "green4:dot:1:green4:square:0.5:green4",
	   "mpfun90" => "green4:dash:1:green4:square:0.5:none",
	   "nag-cmplx" => "red:dot:1:red:square:0.5:red",
	   "nag-workspace" => "red:dash:1:red:square:0.5:none",
	   "nag-inplace" => "red:dot:1:red:triangle-right:0.5:red",
	   "nag-multiple" => "red:dash:1:red:circle:0.5:none",
	   "napack" => "green4:dot:2:green4:none:0.5:none",
	   "nr-c" => "brown:dot:1:brown:square:0.5:none",
	   "nr-f" => "brown:dotdotdash:1:brown:square:0.5:brown",
	   "numutils" => "cyan:dash:1:cyan:triangle-down:0.5:none",
	   "ooura-sg" => "magenta:solid:1:magenta:square:0.7:none",
	   "ooura-sgf" => "magenta:solid:1:magenta:square:0.7:magenta",
	   "ooura-sg2d" => "magenta:solid:1:magenta:square:0.7:none",
	   "ooura-sg2df" => "magenta:solid:1:magenta:square:0.7:magenta",
	   "ooura-sg3d" => "magenta:solid:1:magenta:square:0.7:none",
	   "ooura-sg3df" => "magenta:solid:1:magenta:square:0.7:magenta",
	   "ooura-8g" => "magenta:dot:1:magenta:circle:0.5:none",
	   "ooura-8gf" => "magenta:dot:1:magenta:circle:0.5:magenta",
	   "ooura-4g" => "magenta:dash:1:magenta:triangle-up:0.5:none",
	   "ooura-4gf" => "magenta:dash:1:magenta:triangle-up:0.5:magenta",
	   "ooura-4f2d" => "magenta:dash:1:magenta:triangle-up:0.5:none",
	   "ooura-4f2df" => "magenta:dash:1:magenta:triangle-up:0.5:magenta",
	   "qft" => "violet:solid:1:violet:plus:0.5:none",
	   "ransom" => "green4:solid:1:green4:diamond:0.5:none",
	   "rmayer-ctrig" => "orange:solid:1:orange:circle:0.25:none",
	   "rmayer-buneman" => "orange:dot:1:orange:circle:0.5:none",
	   "rmayer-buneman2" => "orange:dot:1:orange:triangle-up:0.5:none",
	   "rmayer-buneman3" => "orange:dot:1:orange:star:0.5:none",
	   "rmayer-simple" => "orange:dash:1:orange:diamond:0.5:orange",
	   "rmayer-unstable" => "orange:dash:2:orange:diamond:0.5:none",
	   "rmayer-lookup" => "orange:solid:1:orange:triangle-down:0.5:none",
	   "sciport" => "turquoise:solid:1:turquoise:triangle-down:0.4:none",
	   "sgimath" => "black:solid:1:black:star:0.5:none",
	   "singleton" => "red:solid:1:red:triangle-up:0.5:none",
	   "singleton-3d" => "red:solid:1:red:triangle-up:0.5:none",
	   "sorensen-ctfftsr" => "indigo:dot:2:indigo:none:0.25:none",
	   "sorensen-rsrfft" => "indigo:solid:1:indigo:diamond:0.5:indigo",
	   "sorensen-sfftfu" => "indigo:dash:1:indigo:x:0.5:none",
	   "spiral-fft" => "maroon:dot:1:maroon:star:0.5:none",
	   "sunperf" => "black:solid:1:black:star:0.5:none",
	   "temperton" => "grey:solid:2:black:circle:0.25:none",
	   "teneyck" => "grey:dash:2:black:square:0.5:grey",
	   "valkenburg" => "cyan:solid:1:cyan:star:0.5:none",
	   "vbigdsp" => "red:solid:1:red:star:0.5:none",
	   "vdsp" => "black:solid:1:black:star:0.5:none",
	   "vDSP" => "black:solid:1:black:star:0.5:none",
	   );

%symbols = ("none" => 0, "circle" => 1, "square" => 2, "diamond" => 3,
	    "triangle-up" => 4, "triangle-left" => 5, "triangle-down" => 6,
	    "triangle-right" => 7, "plus" => 8, "x" => 9, "star" => 10, 
	    "char" => 11);

%linestyles = ("none" => 0, "solid" => 1, "dot" => 2, "dash" => 3,
	       "dash2" => 4, "dotdash" => 5, "dotdash2" => 6,
	       "dotdotdash" => 7, "dotdashdash" => 8);

sub printstyle {
    my $style = shift;
    my $setnum = shift;
    my ($linecolor,$linestyle,$linewidth,$symbolcolor,$symbolshape,
	$symbolsize,$symbolfillcolor) = split(/:/,$style);
    if ($linecolor eq "none" || $linestyle eq "none") {
	print "@ s$setnum linestyle ",$linestyles{"none"},"\n";
    }
    else {
	print "@ s$setnum linestyle $linestyles{$linestyle}\n";
	print "@ s$setnum linewidth $linewidth\n";
	print "@ s$setnum line color \"$linecolor\"\n";
    }
    if ($symbolcolor ne "none") {
	# note: '<char>' symbols are not yet implemented here.
	print "@ s$setnum symbol $symbols{$symbolshape}\n"; 
	print "@ s$setnum symbol color \"$symbolcolor\"\n";
	print "@ s$setnum symbol size $symbolsize\n";
	if ($symbolfillcolor ne "none") {
	    print "@ s$setnum symbol fill pattern 1\n";
	    print "@ s$setnum symbol fill color \"$symbolfillcolor\"\n";
	}
    }
}

#############################################################################

# Put the legend in a good(?) place:
print "@ legend char size 0.75\n";
print "@ legend 1.02,0.85\n";
print "@ view xmax 1.0\n"; # make space for the legend

# print "@ xaxis label \"transform size\"\n";
if ($accuracy) {
    print "@ yaxis label \"relative rms error\"\n";
    $best_val = -1;
    $worst_val = 0;
}
else {
    print "@ yaxis label \"speed (mflops)\"\n";
    $best_val = 0;
    $worst_val = 1e20;
}

# Collect the data:
while (<>) {
  if ($accuracy) {
      ($nam, $prob, $siz, $L1f, $L2f, $Linff, $L1b, $L2b, $Linfb) = split / /;
      $val = -$L2f;
  }
  else {
      ($nam, $prob, $siz, $mflops, $tim, $setup_tim) = split / /;
      $val = $mflops;
  }
  next if ($val == 0.0 || $val eq "");
  $tot = $siz;
  $tot =~ s/x/*/g;
  $tot = eval($tot);

  next if ($accuracy and $tot < 8);

  $tots{$siz} = $tot;

  if (! exists($best_vals{$siz}) || $best_vals{$siz} < $val) {
      $best_vals{$siz} = $val;
  }
  if (! exists($worst_vals{$siz}) || $worst_vals{$siz} > $val) {
      $worst_vals{$siz} = $val;
  }
  if ($val > $best_val) { $best_val = $val };
  if ($val < $worst_val) { $worst_val = $val };

  if (! exists($problems{$nam})) { $problems{$nam} = $prob; }
  elsif ($problems{$nam} ne $prob) { $problems{$nam} = "several"; }

  $transform = "$nam:$prob";
  if (! exists($results{$transform})) { 
      $results{$transform} = "$siz:$val";
  }
  else {
      $results{$transform} = $results{$transform} . " $siz:$val";
  }
}

# Set x axis ticks and labels:
print "@ xaxis ticklabel angle 270\n";
print "@ xaxis ticklabel type spec\n";
print "@ xaxis tick type spec\n";
@sizes = sort { $tots{$a} - $tots{$b} } keys(%tots);
print "@ xaxis tick spec ", 1 + $#sizes, "\n";
$labelsize = 30.0 / (1 + $#sizes);
if ($labelsize < 1.0) { print "@ xaxis ticklabel char size $labelsize\n"; }
$ticknum = 0;
foreach $siz (@sizes) {
    print "@ xaxis tick major $ticknum, $ticknum\n";
    print "@ xaxis ticklabel $ticknum, \"",$siz,"\"\n";
    $xval{$siz} = $ticknum;
    $ticknum = $ticknum + 1;
}

# Find the y axis scale from $best_val:
if ($accuracy) {
    $val_increment = 1;
    $best_val = -log(-$best_val) / log(10) + 1.0;
    $worst_val = -log(-$worst_val) / log(10) - 1.0;
}
else {
    $val_increment = 100; # increment for y-axis labels
    if ($best_val > 1500) { $val_increment = 500; }
}
$best_val =~ s/\..*//;
$best_val = $best_val + $val_increment - 1 - ($best_val + $val_increment - 1) % $val_increment;
$worst_val =~ s/\..*//;
$worst_val = $worst_val + $val_increment - 1 - ($worst_val + $val_increment - 1) % $val_increment;

# Set axis scales, grids, and colors:
print "@ world xmin 0\n";
print "@ world xmax ",$ticknum-1,"\n";
if ($accuracy) {
    $worst_val = $best_val - 2 if ($accurate_only);
    $best_val = exp(-log(10) * $best_val);
    $worst_val = exp(-log(10) * $worst_val);
    print "@ world ymin $best_val\n";
    print "@ world ymax $worst_val\n";
    print "@ yaxes scale Logarithmic\n";
    $val_increment = 10;
}
else {
    print "@ world ymin 0\n";
    print "@ world ymax $best_val\n";
}
print "@ yaxis tick major $val_increment\n";
print "@ yaxis tick minor color \"grey\"\n";
print "@ yaxis tick major color \"grey\"\n";
print "@ yaxis tick major grid on\n";
print "@ xaxis tick major color \"grey\"\n";
print "@ xaxis tick major grid on\n";
print "@ autoscale onread none\n";  # requires recent version of Grace

# add grid lines?

# Make sure results are sorted in increasing order of xval, or grace will
# connect the dots strangely:
foreach $transform (keys %results) {
    $results{$transform} = join(" ", sort { 
	($siz_a, $val_a) = split(/:/,$a);
	($siz_b, $val_b) = split(/:/,$b);
	$xval{$siz_a} - $xval{$siz_b};
    } split(/ /,$results{$transform}));
}

# Compute median "normalized value" of each transform, for sorting purposes:
@norm_vals = ();
foreach $transform (keys %results) {
    @nvs = ();
    foreach $speed (split(/ /, $results{$transform})) {
	($siz, $val) = split(/:/,$speed);
	if ($accuracy) {
	    $val = log(-$val) / log(-$best_vals{$siz});
	}
	else {
	    $val = $val / $best_vals{$siz};
	}
	$nvs[$#nvs + 1] = $val;
    }
    @nvs = sort { 100000 * ($b - $a) } @nvs;
    $norm_val = $nvs[($#nvs + 1) / 2];
    $norm_val = $norm_val - 1e-10 if (exists($transforms{$norm_val}));
    $norm_vals[$#norm_vals + 1] = $norm_val;
    $transforms{$norm_val} = $transform;
}

# Print out the data sets for each transform, in descending order of "value":
$setnum = 0;
%done = ();
foreach $norm_val (sort { 100000 * ($b - $a) } @norm_vals) {
    $transform = $transforms{$norm_val};
    ($nam, $prob) = split(/:/,$transform);
    $namleg = $nam;

    $namleg = "fxt-matrixfft" if ($namleg eq "fxt-4step");

    # get transform "family"
    ($nam0,$namrest) = split(/\-|:|77|90/, $nam);
    $nam0 = $nam if ($nam eq "intel-ipps");
    $nam0 = $nam if ($nam eq "fftw3-r2r");
    $nam0 = $nam if ($namleg eq "fxt-matrixfft");
    $nam0 = $nam if ($nam0 eq "dsp79");
    if ($nam0 eq "fftw3" || $nam0 eq "fftw3-r2r") {
	if ($prob =~ /..i./) {
	    $nam0 = "$nam0 in-place";
	    $namleg = "$namleg in-place";
	} else {
	    $nam0 = "$nam0 out-of-place";
	    $namleg = "$namleg out-of-place";
	}
    }

    # check if we have already output a relation of this transform
    next if ($no_dups && exists($done{$nam0}));
    $done{$nam0} = 1;

    # legend should include the problem if this transform solves more than
    # one problem in this data:
    if (!$no_dups and $problems{$nam} ne $prob) {
	print "@ s$setnum legend \"$transform\"\n";
    }
    else {
	print "@ s$setnum legend \"$namleg\"\n";
    }

    if (exists($styles{$transform})) {
	printstyle($styles{$transform}, $setnum);
    } elsif (exists($styles{$namleg})) {
	printstyle($styles{$namleg}, $setnum);
    } elsif (exists($styles{$nam})) {
	printstyle($styles{$nam}, $setnum);
    }

    print "@ target s$setnum\n";
    foreach $speed (split(/ /, $results{$transform})) {
	($siz, $val) = split(/:/,$speed);
	$val = -$val if ($accuracy);
	print "$xval{$siz} $val\n";
    }
    print "&\n";

    $setnum = $setnum + 1;
}
