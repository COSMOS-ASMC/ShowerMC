#!/usr//bin/perl
require 5.006;
use warnings;
use strict;
use Getopt::Long;

our $VERSION = '2.02';

our $opt_his = 0;
our $opt_debug;
our $opt_dbg;
our $out_filename;
GetOptions(
    'his:s' => \$opt_his,
    'debug' => \$opt_debug,
    'dbg'   => \$opt_dbg,
    'f=s'   => \$out_filename,
);

sub usage {
    print <<HERE;
Syntax:  add.pl -f outputfile inputfiles

further options:  -debug   gives debug high level
                  -dbg     gives debug low level

Version: $VERSION
HERE
    exit();
}

usage() if not @ARGV;

# compile a few regular expressions in advance
# NOTE: We explicitly allow d or D in place of the usual eE for the scientific notation. Strange format, it is!
our $FloatingNumber = qr/(?: [+-]? (?=\d|\.\d) \d* (?:\.\d*)? (?:[EedD] [+-]? \d+)? )/x;
our $YieldLine      = qr/ te?xt .* (?:I|i|yield|si(?:mu)?[ ])=? \s* ($FloatingNumber) /x;

# open all the input handles
our @input_handles;
{
    my @input_files = @ARGV;
    foreach my $input_file (@input_files) {
        open my $fh, '<', $input_file or die "Could not open input file '$input_file': $!";
        push @input_handles, $fh;
    }
}

my $output_handle;
if ( defined $out_filename) {
    open $output_handle,">", $out_filename
      or die "can't open file '$out_filename' for writing: $!";
}

my $in_array_block   = 0;
my $n_histograms     = 0;
my $print_to_stdout  = 0;
my $n_columns        = 0;

# initialize histo-weights with N ones
my @hw = ( (1) x @input_handles );
my $hwww;

my $lineno = 0;
my $histo_name;
# read first file line by line
while (1) {
    my ($preceding_ws, $line) = get_line($input_handles[0]);
    last if not defined $line;

    $lineno++; 

    if ($line =~ /endarray/) {
        # end of histogram: reset histo weight et al
        @hw               = ( (1) x @input_handles );
        $in_array_block   = 0;
        $n_columns        = 0;
        $histo_name       = "";
        $n_histograms++;
        $hwww             = 0;

        # if we have endarray and array on the same line, we assume it's
        # a function, so we skip this line
        if ($line =~ /(?<!get|end)array\b/) {
            # copy to output
            print $output_handle $preceding_ws, $line, "\n";

            # skip ahead all other handles by one line
            foreach my $handle (@input_handles[1..$#input_handles]) {
                my $ignored = <$handle>;
            }
            # skip the rest of the processing and go directly to the
            # next line
            next;
        }
    }
    # remember the current histogram name
    elsif  ($line =~ /name/) {
        $histo_name = $line;
        chomp $histo_name;
    }
    
    # if we're currently within an array block, we do our hist
    if ( $in_array_block ) {
        # This is the OLD logic:
        # y = \sum yi * wi  ;  dy = \sum dyi * wi**1.5 ; wi = hwi/hwww or (yi/dyi)**2
        #
        # NOTE: This block is so long and redundant to unroll all the different edge cases
        # out of the loop as an optimization!

        my $w; # FIXME declared outside the block just for debugging output

        my $x  = (split /\s+/, $line, 2)[0];
        my $y  = 0;
        my $dx = 0;
        my $dy = 0;

        # read one line from each file
        my @lines = ($line, map { (get_line($input_handles[$_]))[1] } 1..$#input_handles);

        my $build_up_hwww = ($hwww==0);

        if ($build_up_hwww) { # all weights were explicitly set zero
            # build up the sum of histoweights as we go
            my $histoweight_sum = 0;
            if ($n_columns == 4) {
                # THIS IS ILLEGAL!
                die "Number of columns set to 4, but histoweight line explicitly set to 0 in all files. This is illegal!";
#                foreach my $lineno (0..$#lines) {
#                    my $this_line = $lines[$lineno];
#                    my ($xi, $yi, $dyi, $wi) = split /\s+/, $this_line, 4;
#                    print "    $xi  $yi  $dyi  $wi" if $print_to_stdout;
#                    $y    += $yi*$wi;
#                    $dy   += ($dyi*$wi)**2;
#                    $histoweight_sum += $wi;
#                }
            }
            elsif ($n_columns == 3) {
                foreach my $lineno (0..$#lines) {
                    my $this_line = $lines[$lineno];
                    my ($xi, $yi, $wi) = split /\s+/, $this_line, 3;
                    print "    $xi  $yi  $wi" if $print_to_stdout;
                    $y               += $yi*$wi;
                    $dy              += $wi; # dy really disguises as wi!!! We want the number of histos there
                    $histoweight_sum += $wi; 
                }
            }
            else { # two cols
                foreach my $lineno (0..$#lines) {
                    my $this_line = $lines[$lineno];
                    my ($xi, $yi) = split /\s+/, $this_line, 2;
                    print "    $xi  $yi" if $print_to_stdout;
                    # wi implied one
                    $y += $yi;
                    $histoweight_sum++;
                }
            }

            # for the weighted mean, divide by the sum of weights
            if ($histoweight_sum > 1.e-12) {
                $y  /= $histoweight_sum; # degrades to unweighted mean in case of two columns
                $dy = sqrt($dy)/$histoweight_sum if $n_columns != 3;
            }
            else {
                $dy = $y = 0;
            }
        } # end hwww is zero
        # at least one histoweight non-zero (explicit or implicit by omission)
        else {
            if ($n_columns == 4) {
                foreach my $lineno (0..$#lines) {
                    my $this_line = $lines[$lineno];
                    my ($xi, $yi, $dxi, $dyi) = split /\s+/, $this_line, 4;
                    print "    $xi  $yi  $dxi $dyi" if $print_to_stdout;
                    my $wi = $hw[$lineno];
                    $y  += $yi*$wi;
                    $dx += ($dxi*$wi)**2;
                    $dy += ($dyi*$wi)**2;
                }
            }
            elsif ($n_columns == 3) {
                foreach my $lineno (0..$#lines) {
                    my $this_line = $lines[$lineno];
                    my ($xi, $yi, $dyi) = split /\s+/, $this_line, 3;
                    print "    $xi  $yi  $dyi" if $print_to_stdout;
                    my $wi = $hw[$lineno];
                    $y  += $yi*$wi;
                    $dy += ($dyi*$wi)**2;
                }
            }
            else { # two cols
                foreach my $lineno (0..$#lines) {
                    my $this_line = $lines[$lineno];
                    my ($xi, $yi) = split /\s+/, $this_line, 2;
                    print "    $xi  $yi" if $print_to_stdout;
                    my $wi = $hw[$lineno];
                    $y  += $yi*$wi;
                }
            }
            # for the weighted mean, divide by the sum of weights
            $y  /= $hwww;
            $dy = sqrt($dy) / $hwww;
            $dx = sqrt($dy) / $hwww;
        } # end hwww not zero

        #.....output of result.......... 
        if ( $n_columns == 4 ) {
            printf $output_handle " %10e %10e %10e %10e \n", $x, $y, $dx, $dy;
            print "--> $x  $y  $dx  $dy \n\n" if $print_to_stdout;
            #print "--> $x  $y  $dx  $dy weight: $w \n\n" if $print_to_stdout;
        } 
        elsif ( $n_columns == 3 ) {
            printf $output_handle " %10e %10e %10e \n", $x, $y, $dy;
            print "--> $x  $y  $dy \n\n" if $print_to_stdout;
            #print "--> $x  $y  $dy weight: $w \n\n" if $print_to_stdout;
        }
        else {
            printf $output_handle " %10e %10e \n", $x, $y;
            print "--> $x  $y  \n\n" if $print_to_stdout;
            #print "--> $x  $y  weight: $w \n\n" if $print_to_stdout;
        }
    } # end if in array block 
    
    # parse a histoweight declaration
    elsif ( $line =~ /histoweight \b \s* ($FloatingNumber)/x ) {
        # set primary histogram weight
        ($hw[0] = $1) =~ tr/[dD]/EE/; # might be 1D+01 instead of 1E+01 (DOH!)
        # set secondary histogram weights
        foreach my $fileno (1 .. @input_handles-1) {
            my ($prec_ws, $secondary_line) = get_line($input_handles[$fileno]);
            # if we have a histoweight line everywhere, set it
            if ($secondary_line =~ /histoweight \b \s* ($FloatingNumber)/x) {
                ($hw[$fileno] = $1) =~ tr/[dD]/dE/; # might be 1D+01 instead of 1E+01 (DOH!)
            }
            # otherwise rewind the line we just read and set to one
            else {
                undo_getline($input_handles[$fileno]);
                $hw[$fileno] = 1;
            }
        }
    }
    # handle the remaining text lines including those that have a "yield" number
    else {
        my $yield  = 0;
        my $nyield = 0;
        my $is_yield_line = ($line =~ $YieldLine);
        if ( $is_yield_line ) {
            $yield  = $1;
            $nyield++;
        }

        foreach my $fileno (1 .. @input_handles-1) {
            my ($prec_ws, $inp) = get_line($input_handles[$fileno]);
            # this covers the case that the first file doesn't have a histoweight line,
            # but some of the others do. Suppose we find a histoweight line in some file,
            # we simply set the weight and skip this file handle ahead by one line
            if ($inp =~ /histoweight \b \s* ($FloatingNumber)/x) {
                ($hw[$fileno] = $1) =~ tr/[dD]/dE/; # might be 1D+01 instead of 1E+01 (DOH!)
                ($prec_ws, $inp) = get_line($input_handles[$fileno]);
            }

            if ( $inp =~ $YieldLine ) {
                $yield += $1;
                $nyield++;
            }
            if ($line ne $inp) { # lines differ! Error!
                # different number columns specified:
                if ($line =~ /array/ and $inp =~ /array/ ) { 
                    (undef, my $ncol1) = split /\s+/, $line;
                    (undef, my $ncol2) = split /\s+/, $inp;
                    die "Differences to file number $fileno in ncol! \n--->$line<--\n--->$inp<--\n"
                      if $ncol1 != $ncol2;
                }
                # die if lines differ that don't (both!) match text/txt/xrange/yrange/array
                elsif ( not( $line =~ /(te?xt|[xy]range|array)/ and $inp =~ /\Q$1\E/ ) ) {
                    die "Differences to file number $fileno in line $lineno! \n--->$line--\n--->$inp<--\n";
                }
            }
        }

        if ( $is_yield_line ) {
            $yield /= $nyield;
            $yield = sprintf("%.6g", $yield);
            #if    ( $yield < 1 )   { $yield=int($yield*10000)/10000; }
            #elsif ( $yield < 100 ) { $yield=int($yield*100)/100; }
            #else                   { $yield=int($yield*10)/10; }
            # replace yield:
            $line =~ s/(te?xt.*(?:[Ii]|yield|si(?:mu)?)=?\s*)$FloatingNumber(?![eE\.\-+0-9])/$1$yield/x
              or warn "Could not replace yield number in output!";
        }
        # start of the array block! (next line will go into the $in_array_block mode)
        elsif ( $line =~ /(?<!get|end)array\s*(\d+)/ ) { # a pure 'array' not preceded by get or end
            $n_columns       = $1;
            $print_to_stdout = 0;
            $in_array_block  = 1;
            
            $hwww = 0;
            $hwww += $_ foreach @hw; # sum of weights

            print $output_handle "histoweight $hwww\n";

            print "$histo_name \n" if $opt_dbg;
            $print_to_stdout = 1
              if $opt_debug and ( $n_histograms == $opt_his or $opt_his == 0 );
            print "adding histo $n_histograms : $histo_name\n"  if $print_to_stdout;
            print "number of columns $n_columns \n" if $print_to_stdout;
        } 
        
        # copy the line to the output file
        print $output_handle $preceding_ws, $line, "\n";
    } # end handling arbitrary other text lines


} # end foreach line in primary file

exit(0);


# read another line from the file handle and
# remove any whitespace from the beginning and end of the line
{
    my $seek_back = 0;
    sub get_line {
        my $fh = shift;
        $seek_back = tell($fh);
        my $line;
        if (defined($line = readline($fh))) {
            my $preceding = '';
            # strip extra white space
            $preceding = $1 if $line =~ s/^(\s+)//;
            $line =~ s/\s+$//;
            return($preceding, $line);
        }
        return();
    }

    sub undo_getline {
        my $fh = shift;
        seek($fh, $seek_back, 0);
        return();
    }
}

