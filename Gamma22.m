function Gamma_22=Gamma22(dvdx,u,v,beta1,rho,A_0,K_r)
kappa=((4^4.5)*K_r*(beta1^2))/((rho^2)*(A_0^2));
Gamma_22=(5/8)*dvdx + kappa*(5*u + 3*v)/((u-v)^5);
end 