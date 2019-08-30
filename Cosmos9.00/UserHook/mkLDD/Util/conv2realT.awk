# you must give mu, fai,  cosz
BEGIN{c=2.998e8; torad=3.141592/180. }
{r=$1; sint =sqrt(1.-cosz*cosz) ;
    if(cosz > 0.999) sint=sqrt(1.-0.999**2);
    t0= r*mu*1.e9/c;
    dt = t0* (cos(fai*torad)+1.)*sint;
    dt0 = t0*cos(fai*torad)*sint;
    print r, $2-dt, dt, $2-dt+dt0
}



