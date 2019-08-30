function g2mb(A,x){
    return A/6.022e-4/x
}
NR>1 {sinelp=g2mb($1,$2); stp=g2mb($1, $3); \
    sinelc=g2mb($1,$4); stc=g2mb($1, $5); \
    print $1,$2, $3, sinelp, stp, $4,$5, sinelc, stc}
