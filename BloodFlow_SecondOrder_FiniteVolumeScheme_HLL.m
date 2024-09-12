%clear all
%clc
L=0.06;              % length of the artery (with stenosis at x=L)
n=50;               % number of elements or (equi-lenght) intervals within [0 L] 
m=1000050;
dt=0.00001;        % time step (or delta t)
dx=L/n;             %space discretization length or delta x

Wo=2*pi;
K_r=8*pi*0.0035;    % friction coefficient term (blood viscosity taken is 0.0035 Pa sec or 0.035 poise)
rho=1060;           % blood density 
R_0=0.0055;         % Reference radius R_0
R_L=0.0055;         % Radius of the curvature  and x=L
h=0.0005;           % thickness of the artery
b=4/3;
E=4*100000;           % Young's modulus
R_T=100000000;        % Parameter R_T (total terminal resistance [Steriopulous]
K_s=1.52;             % Parameter K_s 
A_0=pi*(R_0^2);       % Reference Area A_0=pi*R_0^2 
A_s=pi*(R_L^2);       % Area A_s=pi*R_L^2
beta1=h*E*b*sqrt(pi); % Beta =h*E*b*sqrt(pi)
epsln=2e-5 + 0.86393e-5;     
x1=linspace(0, L, n);
Qbound0=zeros(m,1);
%% Harmonics of sine and cosine for the fourier series of the inlet $Q_{in} 

%% Initialize cross-section area at x=0 for U_1 
Abound01=zeros(m,1);
%Abound02=zeros(m,1);

%% Initialize cross-section area at x=L for U_1 
AboundL1=zeros(m,1);
%AboundL2=zeros(m,1);

x=linspace(0, L, n+1);            % the lenght L (in meters) is divided into n components,
gamma=4*sqrt(beta1/(2*rho*A_0));

%% Initialize vectors U, u and Q_predicted for first initial condition [U0_11 U0_21]
U_1= zeros(m,2*n);
u_1= zeros(m,n);
Q_predicted1= zeros(m-1,n);

Pressure_drop=zeros(m,1);
V_1_L=zeros(m,1);
%% Initialize $A_{1/2+}$, $A_{1/2-}$, $Q_{1/2+}$, $Q_{1/2-}$ and $u_{1/2+}=Q_{1/2+}/A_{1/2+}$, $u_{1/2+}=Q_{1/2-}/A_{1/2-}$
A1_half_plus=zeros(m-1,n);
A1_half_minus=zeros(m-1,n);

u1_half_plus=zeros(m-1,n);
u1_half_minus=zeros(m-1,n);

Q1_half_plus=zeros(m-1,n);
Q1_half_minus=zeros(m-1,n);
%% Initialize Ghost cell values for at $x=0$ and at $x=L$
Ghost_A1_0=zeros(m-1,1);
Ghost_Q1_0=zeros(m-1,1);
Ghost_u1_0=zeros(m-1,1);

Ghost_A1_L=zeros(m-1,1);
Ghost_Q1_L=zeros(m-1,1);
Ghost_u1_L=zeros(m-1,1);
%% Initialize spatial derivatives of $[A Q]$
dQdx=zeros(m-1,n);
dAdx=zeros(m-1,n);
%% Initialize slope limiter variables for Van Leer slope limiter
Phi_A1=zeros(m-1,n);
Phi_Q1=zeros(m-1,n);
r_A1=zeros(m-1,n);
r_Q1=zeros(m-1,n);
%% Initializing the HLL Flux vectors for U_1 
 FU_11=zeros(m-1,n+1);
 FU_21=zeros(m-1,n+1);

%% Initialize Q at x=0 and x=L for U_1
QboundL1=zeros(m,1);          % intialize Q_in(t)
uboundL1=zeros(m,1);          % initialize velocity vector Q/A 

%% Initialize Q at x=0 and x=L for U_2. Note that Qbound0 is same for both U_1 and U_2 
%QboundL2=zeros(m,1);
%uboundL2=zeros(m,1);          % initialize velocity vector Q/A 

%% Initialize Q at x=0 and x=L for initial time instance t=0 
% (Note that Q_{in}(t) is same for U_1 and U_2) hence same Qbound0 for both
% U_1 and U_2
 % U0_11=Initial_A_new;  
 % U0_21=Initial_Q_new;
 U0_11=Initial_A2_new;  
 U0_21=Initial_Q2_new;
%% Initial condition of vector U=[A Q]' and velocity vector u
U_1(1,:)= [U0_11 U0_21];
u_1(1,:)= U0_21./U0_11;

%% Initialize c_1 and c_2 vectors for U_1 needed for calculating HLL Flux at interfaces of finite volume elements
% For U_1
c_11=zeros(m-1,n-1);
c_21=zeros(m-1,n-1);

c_11G0=zeros(m-1,1);
c_21G0=zeros(m-1,1);
c_11GL=zeros(m-1,1);
c_21GL=zeros(m-1,1);
%% Coefficient for Fourier series needed for $Q_{in}(t)$
B=0.1*[0.13368e-3 -0.12280e-3 0.22459e-4 0.22693e-4 0.22398e-5 -0.22315e-4  0.10065e-4 -0.21066e-5 0.90633e-5  -0.85422e-5];
A=0.1*[-0.88455e-4 -0.52515e-4 0.86471e-4 -0.26395e-4 -0.12987e-4 0.20133e-5 0.70896e-5 0.32577e-5 -0.56573e-5 -0.193023e-5]; 
%% 0th Harmonics of Fourier series
f_0=2e-5 + 0.86393e-5; 
%% $Q_in(dt)$ or inflow at $dt$ time
Qbound0(1)= f_0 + A(1)*cos(Wo*1*dt) + B(1)*sin(Wo*1*dt)+A(2)*cos(Wo*2*dt) + B(2)*sin(Wo*2*dt) + A(3)*cos(Wo*3*dt) + B(3)*sin(Wo*3*dt)+ A(4)*cos(Wo*4*dt) + B(4)*sin(Wo*4*dt) + A(5)*cos(Wo*5*dt) + B(5)*sin(Wo*5*dt)... 
 + A(6)*cos(Wo*6*dt) + B(6)*sin(Wo*6*dt) + A(7)*cos(Wo*7*dt) + B(7)*sin(Wo*7*dt) + A(8)*cos(Wo*8*dt) + B(8)*sin(Wo*8*dt) + A(9)*cos(Wo*9*dt) + B(9)*sin(Wo*9*dt) + A(10)*cos(Wo*10*dt) + B(10)*sin(Wo*10*dt);  

%% Variables calculated recursively for all time t=m*dt
for j=2:m
 t=dt*j;
 %% Update Boundary condition Q(0,t)
  Qbound0(j)= f_0 + A(1)*cos(Wo*1*dt*j) + B(1)*sin(Wo*1*dt*j)+A(2)*cos(Wo*2*dt*j) + B(2)*sin(Wo*2*dt*j) + A(3)*cos(Wo*3*dt*j) + B(3)*sin(Wo*3*dt*j)+ A(4)*cos(Wo*4*dt*j) + B(4)*sin(Wo*4*dt*j) + A(5)*cos(Wo*5*dt*j) + B(5)*sin(Wo*5*dt*j)... 
 + A(6)*cos(Wo*6*dt*j) + B(6)*sin(Wo*6*dt*j) + A(7)*cos(Wo*7*dt*j) + B(7)*sin(Wo*7*dt*j) + A(8)*cos(Wo*8*dt*j) + B(8)*sin(Wo*8*dt*j) + A(9)*cos(Wo*9*dt*j) + B(9)*sin(Wo*9*dt*j) + A(10)*cos(Wo*10*dt*j) + B(10)*sin(Wo*10*dt*j);  

 %% Calculating variables for previous time instant 
 %% Since at x=0 Q(0) is known we use the following scheme

 %% Specify cnear and unear for U_1
  cnear01=c_variable(U_1(j-1,1), rho, beta1, A_0); 
  unear01=(U_1(j-1,n+1)/U_1(j-1,1));

  %% Solve the polynomial to get A(0)

%% Caluclate Abound0 for U_1
  polynomial_01=[gamma^4 -(4*cnear01-unear01)^4 -4*Qbound0(j-1)*((4*cnear01-unear01)^3) -6*(Qbound0(j-1)^2)*((4*cnear01-unear01)^2) 4*(Qbound0(j-1)^3)*(unear01 - 4*cnear01) -Qbound0(j-1)^4];
  
  root_01=roots(polynomial_01);
  root_01=root_01(imag(root_01)==0);
  dist_01=zeros(1,length(root_01));
   if length(root_01)>1

       for q_01=1:length(root_01)

             dist_01(q_01)=abs(root_01(q_01)- U_1(j-1,1));
       end
       
         f01= find(dist_01==min(dist_01));
         Abound01(j-1)=root_01(f01); 
      
    else

      Abound01(j-1)=root_01;
   end
     Ghost_A1_0(j-1)=Abound01(j-1);
     Ghost_Q1_0(j-1)=Qbound0(j-1);
     Ghost_u1_0(j-1)=Ghost_Q1_0(j-1)/Ghost_A1_0(j-1);

  %% Since at x=L Q(L) (as it is a conserved quantity) is known we use the following scheme
   % cnearL and unearL for U_1
   cnearL1=c_variable(U_1(j-1,n),rho,beta1,A_0);
   unearL1=(U_1(j-1,2*n)/(U_1(j-1,n)));

  %% Solve the polynomial to get A(L) for U_1
 % Declare coefficients of the polynomial bar_d1*A^(5/2)=bar_a1A^2 + bar_b1*A +bar_c1, for the BC a x=L
  bar_a1=(beta1/sqrt(A_0)) + R_T*U_1(j-1,2*n) + (0.5*K_s*rho*U_1(j-1,2*n)^2)/(A_s^2);
  bar_b1=-(K_s*rho*U_1(j-1,2*n)^2)/A_s;
  bar_c1=(K_s*rho*0.5)*(U_1(j-1,2*n)^2);
  bar_d1=(beta1/A_0);
  
  
  polynomial_L1=[-bar_d1^2 bar_a1^2 2*bar_a1*bar_b1 (bar_b1^2+(2*bar_c1*bar_a1)) 2*bar_b1*bar_c1 (bar_c1^2)];
  root_L1=roots(polynomial_L1);
  root_L1=root_L1(imag(root_L1)==0);
  dist_L1=zeros(1,length(root_L1));

  if length(root_L1)>1

       for q_L1=1:length(root_L1)
             dist_L1(q_L1)=abs(root_L1(q_L1)- U_1(j-1,n));
       end
         fL1= find(dist_L1==min(dist_L1));
         AboundL1(j-1)=root_L1(fL1);  


  else

      AboundL1(j-1)=root_L1;
  end
    Ghost_A1_L(j-1)=AboundL1(j-1);
    Ghost_Q1_L(j-1)=U_1(j-1,2*n);
    Ghost_u1_L(j-1)=Ghost_Q1_L(j-1)/Ghost_A1_L(j-1);

   

  %% Define the flux at the boundary edges  for U_1
   %% Flux for the 1st interface (i.e., 1-1/2) based on Ghost cell and first cell
  %% Flux at x=0
  % For U_1
  c_11G0(j-1)= min([lambda1(Ghost_Q1_0(j-1), Ghost_A1_0(j-1),rho,beta1,A_0), lambda1(U_1(j-1,n+1),U_1(j-1,1),rho,beta1,A_0), lambda2(Ghost_Q1_0(j-1), Ghost_A1_0(j-1),rho,beta1,A_0), lambda2(U_1(j-1,n+1),U_1(j-1,1),rho,beta1,A_0)]);
  c_21G0(j-1)= max([lambda1(Ghost_Q1_0(j-1), Ghost_A1_0(j-1),rho,beta1,A_0), lambda1(U_1(j-1,n+1),U_1(j-1,1),rho,beta1,A_0), lambda2(Ghost_Q1_0(j-1), Ghost_A1_0(j-1),rho,beta1,A_0), lambda2(U_1(j-1,n+1),U_1(j-1,1),rho,beta1,A_0)]);
  
  FU_11(j-1,1)= (c_21G0(j-1)*Ghost_Q1_0(j-1) - c_11G0(j-1)*U_1(j-1,n+1))/(c_21G0(j-1)-c_11G0(j-1)) + ((c_11G0(j-1)*c_21G0(j-1))/(c_21G0(j-1)-c_11G0(j-1)))*(U_1(j-1,1) - Ghost_A1_0(j-1));
  FU_21(j-1,1)= (c_21G0(j-1)*Flux(Ghost_Q1_0(j-1), Ghost_A1_0(j-1), rho, beta1, A_0) - c_11G0(j-1)*Flux(U_1(j-1,n+1), U_1(j-1,1), rho, beta1, A_0))/(c_21G0(j-1)-c_11G0(j-1)) + ((c_11G0(j-1)*c_21G0(j-1))/(c_21G0(j-1)-c_11G0(j-1)))*(U_1(j-1,n+1) - Ghost_Q1_0(j-1));
  %% Spatial Derivative of A for first cell
  dAdx(j-1,1)=(U_1(j-1,1)-Ghost_A1_0(j-1))/dx;
  
  %% Second-order reconstruction without minmod 
  %% First we define slope limiters (van Leer)
  r_A1(j-1,1)=(U_1(j-1,1)-Ghost_A1_0(j-1))/(U_1(j-1,2)-U_1(j-1,1));
  r_Q1(j-1,1)=(U_1(j-1,n+1)-Ghost_Q1_0(j-1))/(U_1(j-1,n+2)-U_1(j-1,n+1));
  
  Phi_A1(j-1,1)=(r_A1(j-1,1)+abs(r_A1(j-1,1)))/(1+abs(r_A1(j-1,1)));
  Phi_Q1(j-1,1)=(r_Q1(j-1,1)+abs(r_Q1(j-1,1)))/(1+abs(r_Q1(j-1,1)));
  
  %% Next we define the interface values of state variables
  A1_half_plus(j-1,1)=U_1(j-1,1) + (1/2)*Phi_A1(j-1,1)*(U_1(j-1,2)-U_1(j-1,1));
  A1_half_minus(j-1,1)=U_1(j-1,1) - (1/2)*Phi_A1(j-1,1)*(U_1(j-1,2)-U_1(j-1,1));

  Q1_half_plus(j-1,1)=U_1(j-1,n+1) + (1/2)*Phi_Q1(j-1,1)*(U_1(j-1,n+2)-U_1(j-1,n+1));
  Q1_half_minus(j-1,1)=U_1(j-1,n+1) - (1/2)*Phi_Q1(j-1,1)*(U_1(j-1,n+2)-U_1(j-1,n+1));

  
 
  %% Flux at x=L
   c_11GL(j-1)= min([lambda1(Ghost_Q1_L(j-1), Ghost_A1_L(j-1),rho,beta1,A_0), lambda1(U_1(j-1,2*n),U_1(j-1,n),rho,beta1,A_0), lambda2(Ghost_Q1_L(j-1), Ghost_A1_L(j-1),rho,beta1,A_0), lambda2(U_1(j-1,2*n),U_1(j-1,n),rho,beta1,A_0)]);
   c_21GL(j-1)= max([lambda1(Ghost_Q1_L(j-1), Ghost_A1_L(j-1),rho,beta1,A_0), lambda1(U_1(j-1,2*n),U_1(j-1,n),rho,beta1,A_0), lambda2(Ghost_Q1_L(j-1), Ghost_A1_L(j-1),rho,beta1,A_0), lambda2(U_1(j-1,2*n),U_1(j-1,n),rho,beta1,A_0)]);
 
   FU_11(j-1,n+1)= (c_21GL(j-1)*U_1(j-1,2*n) - c_11GL(j-1)*Ghost_Q1_L(j-1))/(c_21GL(j-1)-c_11GL(j-1)) + ((c_11GL(j-1)*c_21GL(j-1))/(c_21GL(j-1)-c_11GL(j-1)))*(Ghost_A1_L(j-1) - U_1(j-1,n));
   FU_21(j-1,n+1)= (c_21GL(j-1)*Flux(U_1(j-1,2*n), U_1(j-1,n), rho, beta1, A_0) - c_11GL(j-1)*Flux(Ghost_Q1_L(j-1), Ghost_A1_L(j-1), rho, beta1, A_0))/(c_21GL(j-1)-c_11GL(j-1)) + ((c_11GL(j-1)*c_21GL(j-1))/(c_21GL(j-1)-c_11GL(j-1)))*(Ghost_Q1_L(j-1) - U_1(j-1,2*n));
   dAdx(j-1,n)=(Ghost_A1_L(j-1)-U_1(j-1,n))/dx;
 
   
  %% First we define slope limiters (van Leer)
   r_A1(j-1,n)=(U_1(j-1,n)-U_1(j-1,n-1))/(Ghost_A1_L(j-1)-U_1(j-1,n));
   Phi_A1(j-1,n)=(r_A1(j-1,n)+abs(r_A1(j-1,n)))/(1+abs(r_A1(j-1,n)));
   Phi_Q1(j-1,n)=2;

  %% Next we define the interface values of state variables
  A1_half_plus(j-1,n)=U_1(j-1,n) + (1/2)*Phi_A1(j-1,n)*(Ghost_A1_L(j-1)-U_1(j-1,n));
  A1_half_minus(j-1,n)=U_1(j-1,n) - (1/2)*Phi_A1(j-1,n)*(Ghost_A1_L(j-1)-U_1(j-1,n));

  Q1_half_plus(j-1,n)=U_1(j-1,2*n) + (1/2)*Phi_Q1(j-1,n)*(Ghost_Q1_L(j-1)-U_1(j-1,2*n));
  Q1_half_minus(j-1,n)=U_1(j-1,2*n) - (1/2)*Phi_Q1(j-1,n)*(Ghost_Q1_L(j-1)-U_1(j-1,2*n));

  QboundL1(j-1)=U_1(j-1,2*n);
  
 %% Assign values of variables for each interval i in the space x.   
 for i=2:n
   
%% Define the values c1 and c2 for U_1

   if i<n
    %% Second-order reconstruction without minmod
    %% First we define slope limiters (van Leer)
    r_A1(j-1,i)=(U_1(j-1,i)-U_1(j-1,i-1))/(U_1(j-1,i+1)-U_1(j-1,i));
    r_Q1(j-1,i)=(U_1(j-1,n+i)-U_1(j-1,n+i-1))/(U_1(j-1,n+i+1)-U_1(j-1,n+i));
    
    Phi_A1(j-1,i)=(r_A1(j-1,i)+abs(r_A1(j-1,i)))/(1+abs(r_A1(j-1,i)));
    Phi_Q1(j-1,i)=(r_Q1(j-1,i)+abs(r_Q1(j-1,i)))/(1+abs(r_Q1(j-1,i)));

    %% Next we define the interface values of state variables
     A1_half_plus(j-1,i)=U_1(j-1,i) + (1/2)*Phi_A1(j-1,i)*(U_1(j-1,i+1)-U_1(j-1,i));
     A1_half_minus(j-1,i)=U_1(j-1,i) - (1/2)*Phi_A1(j-1,i)*(U_1(j-1,i+1)-U_1(j-1,i));
     
     Q1_half_plus(j-1,i)=U_1(j-1,n+i) + (1/2)*Phi_Q1(j-1,i)*(U_1(j-1,n+i+1)-U_1(j-1,n+i));
     Q1_half_minus(j-1,i)=U_1(j-1,n+i) - (1/2)*Phi_Q1(j-1,i)*(U_1(j-1,n+i+1)-U_1(j-1,n+i));
 
    
   end
    %% Define HLL Flux scheme for U_1 
   c_11(j-1,i-1)= min([lambda1(Q1_half_plus(j-1,i-1), A1_half_plus(j-1,i-1),rho,beta1,A_0), lambda1(Q1_half_minus(j-1,i),A1_half_minus(j-1,i),rho,beta1,A_0), lambda2(Q1_half_plus(j-1,i-1), A1_half_plus(j-1,i-1),rho,beta1,A_0), lambda2(Q1_half_minus(j-1,i),A1_half_minus(j-1,i),rho,beta1,A_0)]);
   c_21(j-1,i-1)= max([lambda1(Q1_half_plus(j-1,i-1), A1_half_plus(j-1,i-1),rho,beta1,A_0), lambda1(Q1_half_minus(j-1,i),A1_half_minus(j-1,i),rho,beta1,A_0), lambda2(Q1_half_plus(j-1,i-1), A1_half_plus(j-1,i-1),rho,beta1,A_0), lambda2(Q1_half_minus(j-1,i),A1_half_minus(j-1,i),rho,beta1,A_0)]);

  if c_11(j-1,i-1)>=0

       FU_11(j-1,i)= Q1_half_plus(j-1,i-1);
       FU_21(j-1,i)= Flux(Q1_half_plus(j-1,i-1), A1_half_plus(j-1,i-1), rho, beta1, A_0); 

  elseif c_11(j-1,i-1)<0 && c_21(j-1,i-1)>0

       FU_11(j-1,i)= (c_21(j-1,i-1)*Q1_half_plus(j-1,i-1) - c_11(j-1,i-1)*Q1_half_minus(j-1,i))/(c_21(j-1,i-1)-c_11(j-1,i-1)) + ((c_11(j-1,i-1)*c_21(j-1,i-1))/(c_21(j-1,i-1)-c_11(j-1,i-1)))*(A1_half_minus(j-1,i) - A1_half_plus(j-1,i-1));
       FU_21(j-1,i)= (c_21(j-1,i-1)*Flux(Q1_half_plus(j-1,i-1), A1_half_plus(j-1,i-1), rho, beta1, A_0) - c_11(j-1,i-1)*Flux(Q1_half_minus(j-1,i), A1_half_minus(j-1,i), rho, beta1, A_0))/(c_21(j-1,i-1)-c_11(j-1,i-1)) + ((c_11(j-1,i-1)*c_21(j-1,i-1))/(c_21(j-1,i-1)-c_11(j-1,i-1)))*(Q1_half_minus(j-1,i) - Q1_half_plus(j-1,i-1));

  elseif c_21(j-1,i-1)<=0

        FU_11(j-1,i)= Q1_half_minus(j-1,i);
        FU_21(j-1,i)= Flux(Q1_half_minus(j-1,i), A1_half_minus(j-1,i), rho, beta1, A_0); 

  end 

  %% $dQ/dx$ 
  %dQdx(j-1,i)=FU_11(j-1,i)/dx;

  %% $dA/dx$ 
  dAdx(j-1,i)=(U_1(j-1,i)-U_1(j-1,i-1))/(dx);
 
 end  %% Ending of the For Loop (i=2:n) for number of elements

 %% Update the Area A at j*dt time for U_1
  U_1(j,1:n)= U_1(j-1,1:n) - (dt/dx)*(FU_11(j-1,2:n+1) - FU_11(j-1,1:n));
  dQdx(j-1,1:n)=(FU_11(j-1,2:n+1) - FU_11(j-1,1:n))/dx;
  %dAdx(j-1,1:n)=([U_1(j-1,1:n) Ghost_A1_L(j-1)]-[Ghost_A1_0(j-1) U_1(j-1,1:n-1)])/dx;
 %% Semi-implicit treatment of friction term for U_1
  Q_predicted1(j,1:n)= U_1(j-1,n+1:2*n) - (dt/dx)*(FU_21(j-1,2:n+1) - FU_21(j-1,1:n));

  u_1(j,1:n)= Q_predicted1(j,1:n)./(U_1(j,1:n) + dt*K_r);
  %% Update the Flow Q at j*dt time 
  U_1(j,n+1:2*n)= u_1(j,1:n).*U_1(j,1:n);

  %% Calculate Pressure drop
  V_1_L(m)=(U_1(m,2*n)/AboundL1(m));
  Pressure_drop(m)=PressureDrop(V_1_L(m),AboundL1(m),A_s,K_s,rho);



  

end  % Ending the FOR Loop for time instants j=2:m 

T=linspace(0,t,m);
T1=linspace(0,t-(50*dt),m-50);

 %% Plot Discharge Q over time t
figure(1)
plot(T1, Q1_half_plus(1:m-50,2), T1, U_1(1:m-50,n+12), T1, U_1(1:m-50,n+n/2),T1, U_1(1:m-50,n+37), T1, U_1(1:m-50,2*n - 1), T1, Qbound0(1:m-50),'--', 'LineWidth',2); 
%plot(T1, Q_1(1:m-50,n+1), T1, U_1(1:m-50,n+12), T1, U_1(1:m-50,n+n/2),T1, U_1(1:m-50,n+37), T1, U_1(1:m-50,2*n), T1, Qbound0(1:m-50),'--', 'LineWidth',2); 
xlim([0 10]);
%xlim([0 3]);
fontsize(50,"points")
legend('$Q_1(0^+,t)$','$Q_1(L/4,t)$','$Q_1(L/2,t)$','$Q_1(3L/4,t)$','$Q_1(L,t)$','$Q_{in}(t)$','Interpreter','latex')
xlabel('$t [sec]$','Interpreter','latex');
ylabel('$Q [m^3~sec^{-1}]$','Interpreter','latex');
grid on

%% Plot Area A over time t
figure(2)
plot(T1, U_1(1:m-50,1), T1,  U_1(1:m-50,12), T1,  U_1(1:m-50,n/2), T1,  U_1(1:m-50,37), T1, U_1(1:m-50,n), T1, Abound01(1:m-50),'LineWidth',2); % T, Abound02, '--y', T, Abound01, 'r','LineWidth',2); %x, 0.0000014 + 3.45*0.0000001*sin(w*x),'b');
xlim([0 10]);
%xlim([0 3]);
fontsize(50,"points")
legend('$A_1(0^+,t)$','$A_1(L/4,t)$','$A_1(L/2,t)$','$A_1(3L/4,t)$', '$A_1(L,t)$', '$A_1(0,t)$','Interpreter','latex')
xlabel('$t [sec]$','Interpreter','latex');
ylabel('$A [m^2]$','Interpreter','latex');
grid on

