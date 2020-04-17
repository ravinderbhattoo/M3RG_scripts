function E=buck(r,A,rho,C)

E=A*exp(-r/rho)-C./r.^6;