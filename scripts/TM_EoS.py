from sympy import *

#h, T, C_s = symbols('h, T, C_s')
#
#
#dCsdT     = Derivative( C_s ,T )
#
#h         = 2.5*T + sqrt( 2.25*T**2 + 1 )
#C_s       = T * ( 5*h - 8*T )
#C_s      /= 3*h*( h - T )
#C_s       = sqrt(C_s)
#
#dCsdT_Ans = simplify(diff( C_s,T ))
#
#expr = Eq( dCsdT, dCsdT_Ans )
#
#pprint( expr , use_unicode=True )
#
#
#K, T, A = symbols('K T A')
#
#dTdA = Derivative( T, A )
#T    = A/sqrt( 3*A+1 )
#dTdA_Ans = simplify(diff(T,A))
#
#
#expr = Eq( dTdA, dTdA_Ans )
#pprint( expr , use_unicode=True )
#
#K, rho, A = symbols('K rho A')
#
#dAdrho = Derivative(A, rho)
#
#A = K*rho**(2/3)
#
#dAdrho_Ans = diff(A, rho)
#
#expr = Eq(dAdrho, dAdrho_Ans)
#
#pprint( expr, use_unicode=True )


#Gamma, UsSqr, Us, CsSqr, expr, A, h, T = symbols('Gamma UsSqr Us CsSqr expr A h T')
#
#h = 2.5*T+sqrt(2.25*T**2+1)
#
#CsSqr  = T*( 5*h-8*T )
#CsSqr /= 3*h*( h-T )
#
#UsSqr = CsSqr/( 1-CsSqr )
#
#pprint( simplify(UsSqr), use_unicode=True )

U_zeta, U, gamma_zeta, gamma, gamma_s, rho, U, U_s, C_s, h, T = symbols('U_zeta, U, gamma_zeta gamma gamma_s rho U U_s C_s h T')


ODE = Derivative(U, rho)

ODE = Eq(ODE, gamma*U_s/(gamma_s*rho))

pprint( ODE, use_unicode=True )



U = Eq( gamma_s*U_zeta-U_s*gamma_s )

pprint( U, use_unicode=True )
