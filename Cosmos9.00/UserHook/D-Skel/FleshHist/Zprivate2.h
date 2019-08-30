#include "../../Hist/Z90histc.h"
#include "../../Hist/Z90histo.h"
#include "../../Hist/Z90hist1.h"
! # include "../../Hist/Z90hist2.h"
! # include "../../Hist/Z90hist3.h"

!        timing histogram:    
!                   1       2        3        4         nrbin       
!            *    $ |  $    |   $    |   $    |    $      |   $     
!   rbin(i)         1      2         3        4         nrbin       

          type(histogram1):: tspec(4, nrbin, nfai, nsites)
!                at the center.                                     
          type(histogram1):: tspec0(4, nsites)

          type(histogram1):: rspec(4, nfai, nsites)

