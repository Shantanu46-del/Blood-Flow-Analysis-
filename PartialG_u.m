function DG_u=PartialG_u(Q,A,rho,k)


DG_u=((Q^2)/A) + (k/(3*rho*sqrt(pi)))*(A^(3/2));

end