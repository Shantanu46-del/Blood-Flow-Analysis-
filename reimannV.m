function v_reimann=reimannV(A,Q, rho,beta1,A_0)

v_reimann=Q/A - 2*sqrt((2*beta1)/(rho*A_0))*(A^(0.25));
end