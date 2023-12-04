from sympy import *

h, T, K, C_s, rho, Gamma, gamma = symbols('h T K C_s rho Gamma gamma')

################################################

Target = gamma * ( C_s / rho )
DiffTarget = Derivative( Target, rho );

T = K*rho**(Gamma-1)
h = 1 + (Gamma/(Gamma-1))*T
C_s = sqrt(Gamma*T/h)


Ans = diff( gamma*C_s/rho, rho )
expr = Eq( DiffTarget, simplify(Ans) )
pprint( expr  , use_unicode=True )

################################################

f, C_s, U = symbols('f C_s U')
DiffTarget = Derivative( Target, U );
gamma = sqrt(1+U**2)
Ans = diff( gamma*C_s/rho, U )
expr = Eq( DiffTarget, simplify(Ans) )
pprint( expr  , use_unicode=True )

################################################
