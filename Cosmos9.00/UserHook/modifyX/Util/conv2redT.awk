# you need  mu fai cosz
# convert actual arrival time into reduced time
# input time is a series of (r, t) 
#  (r in m.u, t in ns)
#  reduced time = actual time + dt
# output : r, reduced_time, dt reduced_time'
# the last reduced_time  is not used although
# it can be used in stead of 2nd one.

BEGIN{c=2.998e8; torad=3.141592/180. }
{r=$1; if(cosz >0.999) cosz=0.999; 
    sint=sqrt(1.-cosz*cosz) ;
    dt = r*mu*1.e9/c*(cos(fai*torad)+1.)*sint
    dt0 = r*mu*1.e9/c*cos(fai*torad)*sint;
    print r, $2+dt, dt, $2+dt0
}
