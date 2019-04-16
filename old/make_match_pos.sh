#! /bin/bash
#

awk '{ if ($2 ~ "odi_i") print $5, $6 }' ifirst_tol5.dat > tol5_i.pos ;
awk '{ if ($2 ~ "odi_g") print $5, $6 }' ifirst_tol5.dat > tol5_g.pos ;

awk '{ if ($2 ~ "odi_i") print $5, $6 }' ifirst_tol6.dat > tol6_i.pos ;
awk '{ if ($2 ~ "odi_g") print $5, $6 }' ifirst_tol6.dat > tol6_g.pos ;

awk '{ if ($2 ~ "odi_i") print $5, $6 }' ifirst_tol7.dat > tol7_i.pos ;
awk '{ if ($2 ~ "odi_g") print $5, $6 }' ifirst_tol7.dat > tol7_g.pos ;