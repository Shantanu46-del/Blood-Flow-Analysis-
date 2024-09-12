function dGdU=dGdU_cal(u,v,beta1,rho,A_0,K_s,A_s,R_T)
%R_T=100000000;
d1=(R_T*rho^2*A_0^2)/((4^3)*(beta1^2));
d2=(rho^2*A_0^2)/((4^5)*(beta1^2)*A_s);
dGdU=2*rho*(u-v) - 4*d1*((u-v)^3)*(u+v) - d1*((u-v)^4) - 8*K_s*rho*(u+v)*(d2*((u-v)^4)-1)^2 ...
    - 8*K_s*rho*((u+v)^2)*(d2*((u-v)^4)-1)*(4*d2*((u-v)^3)); 
end 