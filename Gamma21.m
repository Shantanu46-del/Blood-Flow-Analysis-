function Gamma_21=Gamma21(dvdx,u,v,beta1,rho,A_0,K_r)
kappa=((4^4.5)*K_r*(beta1^2))/((rho^2)*(A_0^2));
Gamma_21=(3/8)*dvdx - kappa*(3*u + 5*v)/((u-v)^5);
end 