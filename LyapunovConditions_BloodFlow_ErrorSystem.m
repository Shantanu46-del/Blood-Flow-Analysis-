 L=0.06;              % length of the artery (with stenosis at x=L)
 n=50;               % number of elements or (equi-lenght) intervals within [0 L] 
 m=1000050;
 dt=0.00001;         % time step (or delta t)
 t=m*dt;
 mu=4;      
 dx=L/n;             %space discretization length or delta x
 %x=linspace(0, L, n);
 Wo=2*pi;
 K_r=8*pi*0.0035;    % friction coefficient term (blood viscosity taken is 0.0035 Pa sec or 0.035 poise)
 rho=1060;           % blood density 
 R_0=0.0055;         % Reference radius R_0
 R_L=0.0055;         % Radius of the curvature  and x=L
 h=0.0005;           % thickness of the artery
 b=4/3;
 E=4*100000;         %Young's modulus
 R_T=100000000;      % Parameter R_T (total terminal resistance [Steriopulous]
 K_s=1.52;           % Parameter K_s from reference article [Nikos& Seeley]
 A_0=pi*(R_0^2);     % Reference Area A_0=pi*R_0^2 
 A_s=pi*(R_L^2);     % Area A_s=pi*R_L^2
 beta1=h*E*b*sqrt(pi);    %Beta =h*E*b*sqrt(pi)
 %epsln=2e-5 + 0.86393e-5;    
 T1=linspace(0,t-dt*50,m-50);
 x=linspace(0, L, n);
%% Initialize spatail derivates of Riemann variables 
du_reimann=zeros(m-50,n);
dv_reimann=zeros(m-50,n);
 P_1=zeros(1,n);
 P_2=zeros(1,n);
 Dx_P_1=zeros(1,n);
 Dx_P_2=zeros(1,n);
%% Initialize Lambda matrix diagonal elements $\Lambda(x,t)
Lambda_1=zeros(m-50,n);
Lambda_2=zeros(m-50,n);
%% Initialize Gamma matrix elements $\Gamma(x,t)$ 
Gamma_11=zeros(m-50,n);
Gamma_12=zeros(m-50,n);
Gamma_21=zeros(m-50,n);
Gamma_22=zeros(m-50,n);
%% Initalize spatial derivative of $\Lambda(x,t)$
Lambda_1x=zeros(m-50,n);
Lambda_2x=zeros(m-50,n);
%% Initialize Eigenvalues of $R(x,t)$ 
Eig_min_R=zeros(m-50,n);
Eig_max_R=zeros(m-50,n);
min_Eigenvalue_R=zeros(m-50,1);
%% Initialize $b(t)$ and $a(t)$
b_variable=zeros(m-50,1);
a_variable=zeros(m-50,1);
%% Initialize Boundary Lyapunov condition $b(t)^2\frac{|p_2\Lambda_2(L,t)|}{p_1\Lambda_1(L,t}e^{-\mu L}$
b_var_cond=zeros(m-50,1);
%% Initialize $p_1$ and $p_2$ of $\mathcal{P}$ matrix (positive definite matrix in Lyapunov functional)
p1=1;
p2=0.96;
for j=1:m-50
    for i=1:n
     %% Calculating $\mathcal{P(x)}$ matrix  and $\frac{d\mathcal{P}(x)}{dx}$ matrix 
     P_1(i)=p1*exp(-mu*x(i));
     P_2(i)=p2*exp(mu*x(i));
     Dx_P_1(i)=-p1*mu*exp(-mu*x(i));
     Dx_P_2(i)=p2*mu*exp(mu*x(i));
     %% Calculate spatial derivatives of Riemann variables based on $[\frac{d A(x,t)}{dx} \frac{d Q(x,t)}{dx}]^\top$
     du_reimann(j,i)=reimannDU(dAdx(j,i),dQdx(j,i),U_1(j,i),U_1(j,n+i),rho,beta1,A_0);
     dv_reimann(j,i)=reimannDV(dAdx(j,i),dQdx(j,i),U_1(j,i),U_1(j,n+i),rho,beta1,A_0);
     
     %% Calculating Linearized error System Matrix $\Lambda(x,t)$, $\Gamma(x,t)$ and $d(\Lambda(x,t))/dx$
     Lambda_1(j,i)=(5*u_1_reimann(j,i))/8  + (3*v_1_reimann(j,i))/8;
     Lambda_2(j,i)=(3*u_1_reimann(j,i))/8  + (5*v_1_reimann(j,i))/8;

     Lambda_1x(j,i)=(5*du_reimann(j,i))/8  + (3*dv_reimann(j,i))/8;
     Lambda_2x(j,i)=(3*du_reimann(j,i))/8  + (5*dv_reimann(j,i))/8;
  
     Gamma_11(j,i)=Gamma11(du_reimann(j,i),u_1_reimann(j,i),v_1_reimann(j,i),beta1,rho,A_0,K_r);
     Gamma_12(j,i)=Gamma12(du_reimann(j,i),u_1_reimann(j,i),v_1_reimann(j,i),beta1,rho,A_0,K_r);
     Gamma_21(j,i)=Gamma21(dv_reimann(j,i),u_1_reimann(j,i),v_1_reimann(j,i),beta1,rho,A_0,K_r);
     Gamma_22(j,i)=Gamma22(dv_reimann(j,i),u_1_reimann(j,i),v_1_reimann(j,i),beta1,rho,A_0,K_r);
     
     %% Calculating $R(x,t)$ from Lyapunov integral condition corresponding to equation (27) and (36)
     R=([P_1(i) 0;0 P_2(i)]*[Gamma_11(j,i) Gamma_12(j,i);Gamma_21(j,i) Gamma_22(j,i)]) + ([Gamma_11(j,i) Gamma_21(j,i);Gamma_12(j,i) Gamma_22(j,i)]*[P_1(i) 0;0 P_2(i)]) - ([Dx_P_1(i) 0;0 Dx_P_2(i)]*[Lambda_1(j,i) 0;0 Lambda_2(j,i)]) - ([P_1(i) 0;0 P_2(i)]*[Lambda_1x(j,i) 0;0 Lambda_2x(j,i)]);
     Eigenvalue_R=eig(R);
     Eig_min_R(j,i)=min(Eigenvalue_R);
     Eig_max_R(j,i)=max(Eigenvalue_R);
    
    end
    %% Calculating Boundary Lyapunov condition $b(t)^2\frac{|p_2\Lambda_2(L,t)|}{p_1\Lambda_1(L,t}e^{-\mu L}$ corresponding to equation (34)
 b_variable(j)=(dGdU_cal(u_1_reimann(j,n-1),v_1_reimann(j,n),beta1,rho,A_0,K_s,A_s,R_T)/dGdV_cal(u_1_reimann(j,n),v_1_reimann(j,n),beta1,rho,A_0,K_s,A_s,R_T));
 a_variable(j)=-Lambda_2(j,1)/Lambda_1(j,1);

 b_var_cond(j)=-(Lambda_2(j,n)/Lambda_1(j,n))*exp(2*mu*L)*(b_variable(j)^2)*(p2/p1);
 min_Eigenvalue_R(j)=min(Eig_min_R(j));
 
 b_var_cond(j)=-(Lambda_2(j,n)/Lambda_1(j,n))*exp(2*mu*L)*(b_variable(j)^2)*(p2/p1);
end

 % figure(1)
 % plot(T1, b_var_cond, T1, min_Eigenvalue_R, 'LineWidth',3); 
 % fontsize(50,"points")
 % legend('Condition (33)', 'Condition (35)');
 % xlabel('$t [sec]$','Interpreter', 'latex');
 % ylabel('Lyapunov conditions','Interpreter', 'latex');
 % grid on

 figure(2)
 plot(T1, b_var_cond, 'LineWidth',3); 
 fontsize(50,"points")
 %legend('Condition (33)', 'Condition (35)');
 xlabel('$t [sec]$','Interpreter', 'latex');
 ylabel('Lyapunov conditions','Interpreter', 'latex');
 grid on
 % figure(2)
 % plot(T1, a_variable, 'LineWidth',3); 
 % fontsize(50,"points")
 % xlabel('$t [sec]$','Interpreter', 'latex');
 % ylabel('$a(t)$','Interpreter', 'latex');
 % grid on

figure(3)
plot(T1, min_Eigenvalue_R,'LineWidth',3);
fontsize(50,"points")
xlabel('$t [sec]$','Interpreter', 'latex');
ylabel('${\rm min}_{x\in[0, L]}\lambda(R(x,t)$','Interpreter', 'latex');
grid on

 