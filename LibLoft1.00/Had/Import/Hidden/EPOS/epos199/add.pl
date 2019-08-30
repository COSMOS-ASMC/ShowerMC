#!/usr//bin/perl  
require 5.002;
use Getopt::Long;
$opt_his=0;
GetOptions(qw(his:s debug f dbg ));

if ( $#ARGV == -1 ) 
    {
    print "\n";
    print "Syntax:  add.pl -f outputfile inputfiles \n";
    print " \n";
    print " further options:  -debug   gives debug high level\n";
    print "                   -dbg     gives debug low level\n";
    exit;
    } 



if ( $opt_f ) {
    open(out,"> $ARGV[0]") or die "can't open file : $!" ;
    shift;
    $fh=0;
    while ( $#ARGV >= 0 ) {
	for ($j=0;$j<=$#ARGV ;$j++){$backup[$j]=$ARGV[$j]};
	if($#backup > $#ARGV){ $backup[$#backup]=""; $#backup-- } ;
        $name=$ARGV[0];
#	$^I=".bak";              #do a backup of the file
	$^I="";              #do a backup of the file
	$hwgt=0 ; 
	while (<>) {
	    if ( $_ =~ /histoweight/ ){$hwgt=1;}
	    if ( $_ =~ /endarray/ ){$hwgt=0;}
	    elsif ( $_ =~ /array/ && $hwgt == 0 ) {print "histoweight 1. \n" ; }
	    print;
	}
	for ($j=0;$j<=$#backup ;$j++){$ARGV[$j]=$backup[$j]};
        shift;
        if ( -s $name ) { open($fh,$name); print "opening $fh $name\n";  $fh++; } else { print " $name does not exist or is empty \n"; }
        } 
    } 

$flag=0;
$nhis=0;
$line=0;
for($i=0;$i < $fh;$i++) {$hw[$i]=1;}

while (<0>)   #=========================================================================

    {   
    $line++; 
    s/\s+$/\n/;


    if (/endarray/) #-------------------------------------------------------
        {

        for($i=0;$i < $fh;$i++) {$hw[$i]=1;}
        $flag=0;$nhis++;$ncol=0;
        $name="";
        }
    
    
    if  (/name/) #----------------------------------------------------------
        {    
    
        $name=$_;chop $name;
        }
    
    
    if ( $flag ) #----------------------------------------------------------
        {
        
        # y = \sum yi * wi  ;  dy = \sum dyi * wi**1.5 ; wi = hwi/hwww or (yi/dyi)**2  
        #.....first histo............... 
        $col3_zero = 0;
        ($a1,$a2,$a3,$a4)=split;
        print "    $a1  $a2  $a3  $a4" if $pr;
        if($hwww){$w= $hw[$0]/$hwww; if($aa4){$a2 = $a2 *$w; $a3 = $a3 *$w; $a4 = $a4 *$w ; }
		  else{$a2 = $a2 *$w; $a3 = $a3 *$w**1.5;}}
        else{$a2=$a2*$a3  ; }
        #.....other histos...............
        for($i=1;$i < $fh;$i++) 
            {
            $inp = <$i>;
            $inp =~ s/\s+$/\n/;
            ($aa1,$aa2,$aa3,$aa4)=split(' ',$inp); 
            if( $pr ) {print "    $aa1  $aa2  $aa3  $aa4 ";}
            if($hwww){$wi= $hw[$i]/$hwww;
              if($aa4){$a2 += $aa2 *$wi; $a3 += $aa3 *$wi; $a4 += $aa4 *$wi ; }
		      else{$a2 += $aa2 *$wi; $a3 += $aa3 *$wi**1.5;}}
            else{$a2 += $aa2*$aa3 ; $a3 += $aa3 ;}
            }
        #.....output of result.......... 
	if($hwww){  }
        else{if($a3!=0.){$a2 /= $a3 ;}else{$a2=0;}}
         
        if ( $ncol == 4 ) {
            printf out " %10e %10e %10e %10e \n",$a1,$a2,$a3,$a4;
            print "--> $a1  $a2  $a3  $a4 weight: $w \n\n" if $pr  ;
            } 
        elsif ( $ncol == 3 ) {
            printf out " %10e %10e %10e \n",$a1,$a2,$a3;
            print "--> $a1  $a2  $a3 weight: $w \n\n" if $pr  ;
	}
        else {
            printf out " %10e %10e \n",$a1,$a2;
            print "--> $a1  $a2  weight: $w \n\n" if $pr  ;
            }
        } 
    
    
    elsif ( /histoweight/ ) #-----------------------------------------------
        { 
    
        ($dummy,$hw[0],$_)  = split;     
        for($i=1;$i < $fh;$i++) 
            {
            $inp = <$i>;   
            $inp =~ s/\s+$/\n/;
            ($dummy,$hw[$i],$inp) = split(' ',$inp);
            }
            ## for($i=0;$i < $fh;$i++) {print " $i $hw[$i]\n"; }


        }
    else  #------------------------------------------------------------------
        {

        $yield=0;
        $nyield=0;
        if ( /(text|txt).*(I\=|i\=|yield\=|simu |si )\s*([\d\.]+)/ ) {
            $yield= $3;
            $nyield=1;
            }

        for($i=1;$i<$fh;$i++) 
            { 
            $inp=<$i>;  
            $inp =~ s/\s+$/\n/;
            if ( $inp =~ /(text|txt).*(I\=|i\=|yield\=|simu |si )\s*([\d\.]+)/ ) {
                $yield += $3;
                $nyield += 1;
            }
	    if ( $_ ne $inp && ( $_ =~ /array/ && $inp =~ /array/ ) 
	         ){ 
            ($dummy,$ncol1) = split(' ',$inp) ;
            ($dummy,$ncol2) = split(' ',$_) ;
	    if ( $ncol1 ne $ncol2 ){print "differences to file number $i in ncol !!!!!! \n--->$_<--\n--->$inp<--\n";exit}
	}
	    elsif ( $_ ne $inp && !( 
                                  ($_ =~ /text/   && $inp =~ /text/ ) 
                                ||($_ =~ /txt/    && $inp =~ /txt/ )
                                ||($_ =~ /yrange/ && $inp =~ /yrange/ )
                                ||($_ =~ /xrange/ && $inp =~ /xrange/ )
                                ||($_ =~ /array/ && $inp =~ /array/ )
                               ) 
                ){print "differences to file number $i in line $line !!!!!! \n--->$_<--\n--->$inp<--\n";exit }
            }
        if ( /(text|txt).*(I\=|i\=|yield\=|simu |si )\s*([\d\.]+)/ ) {
            if( $yield < 1 ) { $yield=int($yield/$nyield*10000)/10000; }
            elsif( $yield < 100 ) { $yield=int($yield/$nyield*100)/100; }
            else  { $yield=int($yield/$nyield*10)/10; }
            s/(text|txt)(.*)(I\=|i\=|yield\=|simu |si )\s*([\d\.]+)/\1\2\3 $yield/ ;
            }   
        if ( /array/ && ! /endarray/  && ! /getarray/ )
            {
            $flag = 1;
            $hwww=0;
            for($i=0;$i < $fh;$i++) { $hwww+=$hw[$i]; }
            print out "histoweight $hwww\n";
            $pr=0;
            print "$name \n" if $opt_dbg ;
            if( $opt_debug && ( $nhis == $opt_his ||  $opt_his == 0 )){ $pr=1;}
            print "adding histo $nhis : $name\n"  if $pr ;
            ($dummy,$ncol ) = split ;
            print "number of columns $ncol \n" if $pr ;
            } 
        s/ +$//; print out $_;


        }  #---------------------------------------------------------------------


    }









