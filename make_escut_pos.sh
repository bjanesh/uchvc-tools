#! /bin/bash
#

awk '{ if ($11 ~ "odi_i") print $1, $2 }' point_source_info > escut_i.pos 
awk '{ if ($11 ~ "odi_g") print $1, $2 }' point_source_info > escut_g.pos 
