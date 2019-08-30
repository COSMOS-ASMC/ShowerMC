#  this is to get (r, time) data to be used fittng test. r is in Moliere unit
#  time is in nsec.  time is the time at which integral time distribution
#  becomes a given fraction.  Fraction is specified thru nth below.
#  Input is the columns of (rindex, time).
#  Usage: awk -f getTimeFrac.awk nth=..  mu=.. mu0=.. fai=.. cosz=..  file
#         nth=2 --> 5%, nth=3-->10%  nth=4-->20 %
#         mu:  Moliere unit at the place where data is taken
#         mu0: Base Moliere unit where observation is done (say at 875g/cm2)
#         fai: azimuthal angle in deg. of the data
#         cosz: cos of 1ry zenith angle.  
# If don't want to convert time into reduced time make cosz<=0
#  
BEGIN {torad=3.141592/180.; c=2.998e8}
{r = 0.01*10.0**(($1-1)*0.1); t=$nth;
    if(cosz <= 0. ) { dt = 0. }
    else {
	if(cosz > 0.999) {sint = sqrt(1.0-0.999**2)}
	else { sint = sqrt(1.0-cosz**2)}
#                    fai is deg. here
	dt=r*mu*(cos(fai*torad) + 1.0) *sint*1.e9/c;
    }
    scaledr = r*mu/mu0;
    t = t + dt;
    print scaledr, t;
}

