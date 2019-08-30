#!/usr/local/bin/gnuplot -persist
pa(z) = 49.232*z**(-0.05075)+0.7494-0.0188*z
pb(z) = 0.0871*z**(-0.0923)
ra(z) = 124.9*z**(-0.0197)+1.2969
qb(z) =  -0.3459*z**0.05615
a(E,z) = pa(z)*E**(-0.4) + ra(z)
b(E,z) = pb(z)*E**qb(z) + 0.92
xmin(E) = 4*0.511e-3/E
f(x) =x>xmin(E) && x<1 ? 1./x**(2-b(E,z))/(1+a(E,z)*x**b(E,z))**2 : 0
#    EOF

