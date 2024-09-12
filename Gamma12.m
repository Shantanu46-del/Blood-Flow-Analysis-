function Gamma_12=Gamma12(dudx,u,v,beta1,rho,A_0,K_r)
kappa=((4^4.5)*K_r*(beta1^2))/((rho^2)*(A_0^2));
Gamma_12=(3/8)*dudx + kappa*(5*u + 3*v)/((u-v)^5);
end 