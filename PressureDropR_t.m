function pressure=PressureDropR_t(A,Q,beta1,A_0,R_T)
pressure= (beta1/A_0)*(sqrt(A)-sqrt(A_0)) - R_T*Q;
end