!     ****************************************************************
!     *                                                              *
!     * cNeByApproxB:  Electron # and age s by Approox B
!     *                                                              *
! To be updated to be usable upto t = 150 for large E
!     ****************************************************************
!
!           Gives total electron no. and shower age in approximation B
!           Incident particle should be an electron.  2 dim table is used
!

       subroutine cNeByApproxB(jspec, ale, t, eno, s)
       implicit none
       
       integer jspec    ! input.  one of 1, 2, 3
                        !      1: Ne is obtained 
                        !      2, shower age s is obtained
                        !      3, both are obtained
       real*8   ale     ! input.  log10(E0/Ec) where Ec is the critical energy.
                        !         1<ale<11. usabel upto 15 safely.
       real*8    t      ! input.  depth in r%l where Ne is obtained   0<t<40  
                        !         usable upto 55 r.l safely.
       real*8   eno     ! output.  electron no.  in case of jspec=1 or 3
       real*8   s       ! output.  age of shower. in case of jspec=2 or 3
!
!    dddddddddddddddd local local local local ddddddddddddddddddddd
!
!    enot:  containes log10(electron no.) at (t, e) where t runs from
!          0 to 40 step 2 and   e from 1 to 11 step 0.5 ( e is in
!          log10(e0/ec).)
!    stbl: containes  shower age at same (t,e).
!
!  **note** for very small values of eno, they are made to be -9. and
!          corresponding s is 2.5
!
!
      integer  tsize, asize, i
      parameter (tsize = 21,  asize=tsize*tsize)
!
!       ** dimension enot(tsize, tsize), stbl(tsize,tsize) **
      real*8 enot(asize), stbl(asize)
      real*8 ans
!
      data (   enot(i),i=   1,  81)/
     1 0.0  , 0.333, 0.023,-0.411,-0.895,-1.409,-1.940,-2.485,-3.041,
     2-4.573,-6.105,-7.637,-9.000,-9.000,-9.000,-9.000,-9.000,-9.000,
     3-9.000,-9.000,-9.000, 0.0  , 0.741, 0.676, 0.387, 0.004,-0.429,
     4-0.896,-1.385,-1.892,-2.410,-2.942,-3.482,-4.029,-5.561,-7.093,
     5-8.625,-9.000,-9.000,-9.000,-9.000,-9.000, 0.0  , 1.039, 1.190,
     6 1.047, 0.770, 0.418, 0.018,-0.416,-0.873,-1.348,-1.840,-2.345,
     7-2.858,-3.382,-3.912,-4.450,-4.993,-6.525,-8.057,-9.000,-9.000,
     8 0.0  , 1.278, 1.613, 1.610, 1.436, 1.165, 0.832, 0.455, 0.048,
     9-0.385,-0.838,-1.305,-1.787,-2.279,-2.782,-3.292,-3.809,-4.335/
      data (   enot(i),i=  82, 162)/
     1-4.866,-5.401,-5.944, 0.0  , 1.480, 1.975, 2.097, 2.025, 1.835,
     2 1.568, 1.250, 0.890, 0.500, 0.088,-0.345,-0.793,-1.256,-1.729,
     3-2.215,-2.707,-3.210,-3.717,-4.232,-4.754, 0.0  , 1.657, 2.291,
     4 2.529, 2.552, 2.441, 2.240, 1.977, 1.668, 1.321, 0.947, 0.551,
     5 0.134,-0.298,-0.743,-1.203,-1.671,-2.150,-2.635,-3.129,-3.631,
     6 0.0  , 1.817, 2.575, 2.918, 3.030, 2.994, 2.859, 2.651, 2.389,
     7 2.088, 1.752, 1.391, 1.006, 0.605, 0.185,-0.246,-0.690,-1.146,
     8-1.610,-2.082,-2.564, 0.0  , 1.962, 2.832, 3.271, 3.466, 3.503,
     9 3.429, 3.278, 3.065, 2.806, 2.509, 2.182, 1.831, 1.459, 1.068/
      data (   enot(i),i= 163, 243)/
     1 0.662, 0.241,-0.190,-0.633,-1.086,-1.547, 0.0  , 2.095, 3.067,
     2 3.595, 3.868, 3.975, 3.962, 3.863, 3.699, 3.482, 3.225, 2.932,
     3 2.614, 2.270, 1.907, 1.526, 1.131, 0.722, 0.299,-0.132,-0.573,
     4 0.0  , 2.220, 3.286, 3.897, 4.244, 4.414, 4.460, 4.414, 4.295,
     5 4.121, 3.901, 3.645, 3.357, 3.044, 2.707, 2.353, 1.982, 1.594,
     6 1.195, 0.784, 0.360, 0.0  , 2.336, 3.490, 4.178, 4.593, 4.826,
     7 4.929, 4.932, 4.860, 4.726, 4.545, 4.322, 4.067, 3.782, 3.474,
     8 3.144, 2.797, 2.433, 2.054, 1.663, 1.260, 0.0  , 2.446, 3.682,
     9 4.442, 4.922, 5.215, 5.370, 5.424, 5.395, 5.301, 5.157, 4.968/
      data (   enot(i),i= 244, 324)/
     1 4.745, 4.491, 4.208, 3.905, 3.582, 3.239, 2.883, 2.512, 2.127,
     2 0.0  , 2.550, 3.864, 4.691, 5.233, 5.582, 5.789, 5.889, 5.903,
     3 5.851, 5.742, 5.588, 5.394, 5.169, 4.915, 4.636, 4.336, 4.018,
     4 3.680, 3.330, 2.966, 0.0  , 2.649, 4.036, 4.927, 5.528, 5.930,
     5 6.188, 6.332, 6.390, 6.375, 6.303, 6.180, 6.019, 5.821, 5.594,
     6 5.341, 5.064, 4.767, 4.453, 4.121, 3.776, 0.0  , 2.743, 4.200,
     7 5.152, 5.809, 6.262, 6.567, 6.755, 6.853, 6.876, 6.839, 6.751,
     8 6.618, 6.450, 6.248, 6.020, 5.768, 5.492, 5.199, 4.888, 4.561,
     9 0.0  , 2.834, 4.356, 5.367, 6.077, 6.579, 6.929, 7.161, 7.298/
      data (   enot(i),i= 325, 405)/
     1 7.358, 7.354, 7.298, 7.196, 7.055, 6.881, 6.677, 6.447, 6.195,
     2 5.921, 5.631, 5.324, 0.0  , 2.921, 4.506, 5.572, 6.333, 6.883,
     3 7.276, 7.549, 7.724, 7.821, 7.851, 7.826, 7.754, 7.640, 7.492,
     4 7.313, 7.105, 6.875, 6.624, 6.351, 6.063, 0.0  , 3.004, 4.650,
     5 5.770, 6.579, 7.175, 7.610, 7.923, 8.135, 8.266, 8.330, 8.335,
     6 8.292, 8.206, 8.082, 7.928, 7.745, 7.535, 7.304, 7.053, 6.781,
     7 0.0  , 3.085, 4.789, 5.960, 6.815, 7.455, 7.931, 8.282, 8.532,
     8 8.696, 8.792, 8.828, 8.812, 8.753, 8.656, 8.524, 8.364, 8.177,
     9 7.965, 7.733, 7.482, 0.0  , 3.164, 4.922, 6.143, 7.043, 7.725/
      data (   enot(i),i= 406, 441)/
     1 8.242, 8.629, 8.913, 9.111, 9.238, 9.303, 9.315, 9.283, 9.211,
     2 9.103, 8.965, 8.799, 8.609, 8.395, 8.163, 0.0  , 3.239, 5.052,
     3 6.319, 7.264, 7.986, 8.541, 8.964, 9.282, 9.513, 9.670, 9.764,
     4 9.805, 9.797, 9.750, 9.666, 9.549, 9.405, 9.235, 9.042, 8.826/
      data (   stbl(i),i=   1,  81)/
     1 0.0  , 1.199, 1.603, 1.872, 2.074, 2.237, 2.376, 2.497, 2.500,
     2 2.500, 2.500, 2.500, 2.500, 2.500, 2.500, 2.500, 2.500, 2.500,
     3 2.500, 2.500, 2.500, 0.0  , 0.898, 1.267, 1.525, 1.721, 1.878,
     4 2.011, 2.126, 2.227, 2.319, 2.402, 2.479, 2.500, 2.500, 2.500,
     5 2.500, 2.500, 2.500, 2.500, 2.500, 2.500, 0.0  , 0.722, 1.052,
     6 1.296, 1.486, 1.640, 1.770, 1.881, 1.980, 2.069, 2.149, 2.223,
     7 2.291, 2.355, 2.415, 2.471, 2.500, 2.500, 2.500, 2.500, 2.500,
     8 0.0  , 0.608, 0.902, 1.130, 1.313, 1.463, 1.590, 1.700, 1.797,
     9 1.883, 1.962, 2.034, 2.100, 2.162, 2.220, 2.275, 2.326, 2.375/
      data (   stbl(i),i=  82, 162)/
     1 2.422, 2.466, 2.500, 0.0  , 0.527, 0.792, 1.004, 1.178, 1.323,
     2 1.448, 1.556, 1.651, 1.736, 1.814, 1.884, 1.950, 2.010, 2.067,
     3 2.120, 2.171, 2.218, 2.264, 2.307, 2.348, 0.0  , 0.467, 0.707,
     4 0.904, 1.069, 1.210, 1.331, 1.437, 1.531, 1.615, 1.691, 1.761,
     5 1.825, 1.885, 1.941, 1.994, 2.043, 2.090, 2.134, 2.177, 2.217,
     6 0.0  , 0.420, 0.640, 0.823, 0.979, 1.115, 1.232, 1.336, 1.428,
     7 1.512, 1.587, 1.656, 1.720, 1.779, 1.834, 1.886, 1.935, 1.981,
     8 2.025, 2.066, 2.106, 0.0  , 0.383, 0.586, 0.757, 0.905, 1.034,
     9 1.148, 1.249, 1.340, 1.422, 1.497, 1.565, 1.628, 1.687, 1.741/
      data (   stbl(i),i= 163, 243)/
     1 1.792, 1.841, 1.886, 1.930, 1.971, 2.010, 0.0  , 0.353, 0.541,
     2 0.701, 0.841, 0.965, 1.075, 1.174, 1.263, 1.343, 1.417, 1.485,
     3 1.547, 1.605, 1.659, 1.710, 1.758, 1.803, 1.846, 1.887, 1.926,
     4 0.0  , 0.327, 0.503, 0.654, 0.787, 0.905, 1.012, 1.108, 1.194,
     5 1.274, 1.346, 1.413, 1.475, 1.532, 1.586, 1.636, 1.684, 1.728,
     6 1.771, 1.811, 1.850, 0.0  , 0.306, 0.471, 0.613, 0.739, 0.853,
     7 0.955, 1.049, 1.133, 1.211, 1.282, 1.348, 1.409, 1.466, 1.519,
     8 1.569, 1.617, 1.661, 1.703, 1.743, 1.782, 0.0  , 0.287, 0.443,
     9 0.578, 0.698, 0.807, 0.906, 0.996, 1.079, 1.155, 1.225, 1.290/
      data (   stbl(i),i= 244, 324)/
     1 1.350, 1.406, 1.459, 1.509, 1.555, 1.600, 1.642, 1.682, 1.720,
     2 0.0  , 0.271, 0.419, 0.547, 0.661, 0.766, 0.861, 0.949, 1.029,
     3 1.104, 1.172, 1.236, 1.296, 1.352, 1.404, 1.453, 1.499, 1.543,
     4 1.585, 1.625, 1.662, 0.0  , 0.257, 0.397, 0.519, 0.629, 0.729,
     5 0.821, 0.906, 0.984, 1.057, 1.125, 1.187, 1.246, 1.301, 1.353,
     6 1.402, 1.448, 1.491, 1.533, 1.572, 1.610, 0.0  , 0.245, 0.378,
     7 0.495, 0.600, 0.696, 0.785, 0.867, 0.943, 1.014, 1.081, 1.142,
     8 1.200, 1.255, 1.306, 1.354, 1.400, 1.443, 1.484, 1.523, 1.561,
     9 0.0  , 0.233, 0.361, 0.473, 0.573, 0.666, 0.752, 0.832, 0.906/
      data (   stbl(i),i= 325, 405)/
     1 0.975, 1.040, 1.101, 1.158, 1.212, 1.262, 1.310, 1.355, 1.398,
     2 1.439, 1.478, 1.515, 0.0  , 0.223, 0.346, 0.453, 0.550, 0.639,
     3 0.722, 0.799, 0.872, 0.939, 1.003, 1.063, 1.119, 1.172, 1.221,
     4 1.269, 1.314, 1.356, 1.397, 1.435, 1.472, 0.0  , 0.214, 0.332,
     5 0.435, 0.528, 0.614, 0.695, 0.770, 0.840, 0.906, 0.968, 1.027,
     6 1.082, 1.134, 1.183, 1.230, 1.275, 1.317, 1.357, 1.395, 1.432,
     7 0.0  , 0.206, 0.319, 0.418, 0.508, 0.592, 0.669, 0.742, 0.811,
     8 0.875, 0.936, 0.994, 1.048, 1.099, 1.148, 1.194, 1.238, 1.280,
     9 1.320, 1.358, 1.394, 0.0  , 0.198, 0.308, 0.403, 0.490, 0.571/
      data (   stbl(i),i= 406, 441)/
     1 0.646, 0.717, 0.784, 0.847, 0.906, 0.963, 1.016, 1.066, 1.114,
     2 1.160, 1.203, 1.245, 1.284, 1.322, 1.358, 0.0  , 0.191, 0.297,
     3 0.389, 0.473, 0.552, 0.625, 0.694, 0.759, 0.820, 0.878, 0.934,
     4 0.986, 1.036, 1.083, 1.128, 1.171, 1.212, 1.251, 1.289, 1.324/
!
!    dddddddddddddddd local local local local ddddddddddddddddddddd
!
!

 
      if( jspec .gt. 3 .or. jspec .lt. 1)  then
         write(*, *) 'Error input to cNeByApproxB:jspec, ale, t=',
     *                jspec, ale, t
         stop 999
      else
         if(mod(jspec,2) .ne. 0) then
            call k4ptdi(enot, tsize, tsize, tsize, 0.d0, 1.d0, 
     *      2.d0, 0.5d0, t, ale, ans)
            eno=10.0**ans
         endif
         if(mod(jspec/2,2) .ne. 0) then
            call k4ptdi(stbl, tsize, tsize, tsize, 0.d0, 1.d0,
     *      2.d0, 0.5d0, t, ale, ans)
            s=ans
         endif
      endif
      end
