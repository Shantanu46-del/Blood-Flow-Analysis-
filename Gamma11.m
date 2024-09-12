function Gamma_11=Gamma11(dudx,u,v,beta1,rho,A_0,K_r)
kappa=((4^4.5)*K_r*(beta1^2))/((rho^2)*(A_0^2));
Gamma_11=(5/8)*dudx - kappa*(3*u + 5*v)/((u-v)^5);
end 