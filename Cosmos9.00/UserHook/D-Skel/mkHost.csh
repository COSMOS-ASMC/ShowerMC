#!/bin/csh
 @ n = 1
 @ m = $n - 1
 while ( $n < 2049 )
	set numb = ` echo $n | awk '{printf("%4.4d",$1)}' `;
	@ m++;
	echo $m " tasim5"$numb;
	@ m++
	echo $m " tasim5"$numb;
        @ n++;
 end
