function dv_reimann=reimannDV(dAdx, dQdx, A,Q, rho,beta1,A_0)

dv_reimann=(A*dQdx-Q*dAdx)/(A^2) - 0.5*sqrt((2*beta1)/(rho*A_0))*(A^(-0.75))*dAdx;
end