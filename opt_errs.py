import numpy as np

i_sun = 4.58
ap, g, ge, i, ie,  gmi, e_gmi, g_abs, i_abs, mv, MtoL, Ls, m_hi, m_star, His = np.loadtxt('optical_props.txt', usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14), unpack=True)

verr = np.sqrt(ge*ge+ie*ie)
mtol = np.power(10,0.518*gmi-0.152)
mtol_eu = np.power(10,0.518*(gmi+e_gmi)-0.152)
mtol_el = np.power(10,0.518*(gmi-e_gmi)-0.152)
l_star_eu = np.power(10,(i_sun-i_abs-ie)/2.5)
l_star_el = np.power(10,(i_sun-i_abs+ie)/2.5)
m_star_eu = np.log10((l_star_eu)*(mtol_eu))
m_star_el = np.log10((l_star_el)*(mtol_el))
hitostar = m_hi/m_star

for i,a in enumerate(ap):
    print("{:5.1f} {:5.2f} {:4.2f} {:4.2f} {:4.2f} {:4.2f}".format(a,mv[i],verr[i],m_star[i],m_star_eu[i],m_star_el[i]))