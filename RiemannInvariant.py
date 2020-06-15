from sympy import *

h, T, K, Cs, P, A = symbols('h T K Cs P A')


#A = P**(2/3)*K
#Ans = diff(T**(5/3)*(1.5*T+sqrt(2.25*T**2+1)) , T)
#h=2.5*T+sqrt(2.25*T**2+1)
#dh_dT = diff(h,T)
#Ans = integrate( dh_dT/T,T )

Ans = solve( T**(5/3)*(1.5*T+sqrt(2.25*T**2+1)) - A, T)


pprint( simplify( Ans ), use_unicode=True )


#rho     = 1.5*T**2 + T*sqrt( 2.25*T**2 + 1 )
#rho    /= K
#rho     = rho**1.5
#
#h       = 2.5*T + sqrt( 2.25*T**2 + 1 )
#
#Cs      = T * ( 5*h - 8*T )
#Cs     /= 3*h*( h - T )
#Cs      = sqrt(Cs)

#Gamma = symbols('Gamma')
#
#rho = (T/K)**(1/(Gamma-1))
#
#
#drho_dt = diff(rho, T)
#
#h = 1 + (Gamma/(Gamma-1))*T
#
#Cs = Gamma*T/h
#Cs = sqrt(Cs)
#
#integrand = (Cs/rho)*drho_dt
#
#
#pprint( simplify( rho ), use_unicode=True )
#print("======================================")
#pprint( simplify( Cs ), use_unicode=True )
