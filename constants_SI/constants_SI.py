from math import pi

c=299792458            #m/s           EXACT
h=6.62607015-34        #Js=kg m^2/s   EXACT
hbar=h/(2*pi)          #Js=kg m^2/s

e=1.602176634-19       #C             EXACT
m_e=9.10938291e-31     #kg
m_p=1.672621777e-27    #kg

alpha=1/137.035999084  #1
mu_0=2*h*alpha/(c*e*e) #N/A^2
eps_0=1/(mu_0*c*c)     #F/m
epsilon_0=eps_0        #F/m

k_B=1.380649-23        #J/K           EXACT
N_A=6.02214076+23      #particles/mol
N_Avogadro=N_A         #particles/mol

G=6.67384e-11          #N m^2/kg^2 = m^3 /(kg s^2)

### Other useful constants
m_e_eV=m_e*c**2/e      #electron mass in eV
m_p_eV=m_p*c**2/e      #proton   mass in eV
