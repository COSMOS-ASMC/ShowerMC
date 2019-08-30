extern "C" {
#include "Zdef.h"
#include "Zmanagerp.h"
#include "Ztrack.h"
#include "Ztrackp.h"
#include "Ztrackv.h"
#include "Zcode.h"
#include "Zmass.h"
#include "Zprimary.h"
#include "Zprimaryv.h"
#include "Zheavyp.h"
#include "Zincidentv.h"
  /*
c ******************************************************************************
c       The user routines here are used only when you give 161 ~ 199 to the Trace value 
c       so that you  can manage the Cerenkov light output yourself.
c   The main purpose is to enable you to convert each track information to Cerekov light
c   on fly and output it with your desired format.  This will save the output disk file 
c   volume as compared to writing the track infomation directly by the standard way.
c  
c     This file is saved  as ctemplCeren.f
c   Each user hook program is supplied with #include "ctemplCeren.f"
c  (if not supplied, give it somewhere).
c   If you really want to make this file usable, save it as chookCeren.f (or whatever
c   you like), and change its content.
c   Then, you have to change the incldue statement in the user hook program which
c   really uses the chookCeren.f:
c        #include "chookCeren.f"
c ******************************************************************************
c
 */
  void chookcerens_(int *no, primaries *primary, coord *angle){
  /*
    c      implicit none
    c          This is called when one event generation starts.
    c
    c
    c
    integer no  !  input.  Event number.
    record /primaries/ primary  ! input. Primary particle info.
    record /coord/ angle      !  input.  primary angle at the observation depth.
  */

  /*
    c        here  you may put some flag info. as header of each event;
    c     The standard Cerenkov output routine writes the
    c     following:
    c
    c      no      ! event no
    c      primary.particle.code ! intger: partilce code
    c      primary.particle.fm.p(4)  ! energy
    c      angle.r(1), angle.r(2), angle.r(3)      ! direction cos of primary at the observation level.
    c      
    c    
  */
  return ;
}

  void chookceren_(int *ka, int *chrg, double *e1, int *itb,
		   int *it, coord *f, coord *t) {

    return ;
  } 

void chookcerene_(int *ka, int  *chrg, double *e1, int *itb, int *it,
		 coord *f, coord *t){
  return ;
}
}
