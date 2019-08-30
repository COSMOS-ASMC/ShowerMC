# this is to get (r, time) data to be used fittng test. r is in Moliere unit
#  time is in nsec.  time is the time at which integral time distribution
#  becomes a given fraction.  Fraction is specified thru nth below.
#  Input is the columns of (rindex, time).
#  Usage: awk -f getTimeFrac.awk nth=..  mu=.. mu0=.. fai=.. cosz=..  file
#         nth=2 --> 5%, nth=3-->10%  nth=4-->20 %
#         mu:  Moliere unit at the place where data is taken
#         mu0: Base Moliere unit where observation is done (say at 875g/cm2)
#         fai: azimuthal angle in deg. of the data
#         cosz: cos of 1ry zenith angle.  
#    
#  

{r = 0.01*10.0**(($1-1)*0.1); t=$nth;
    if(cosz == 1.0) { dt = 0.}
    else {  sint = sqrt(1.0-cosz**2);
#                    fai is deg. here
	dt=r*mu*(cos(fai*0.001745) + 1.0) *sint*1.e9/2.998e8;
    }
    scaledr = r*mu/mu0;
    t = t + dt;
    print scaledr, t;
}

