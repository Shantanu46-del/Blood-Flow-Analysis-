function u_reimann=reimannU(A,Q, rho,beta1,A_0)

u_reimann=Q/A + 2*sqrt((2*beta1)/(rho*A_0))*(A^(0.25));
end