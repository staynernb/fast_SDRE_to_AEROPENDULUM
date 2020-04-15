%%Teste model_Leis de Newton
%% Aero Pêndulo dual

Jteta = 0.005425;    %1.0348 ;  (1/12*M*((b^2) + c^2)) + M*L1^2 
Jphi = 0.00201;  %0.0451;  (1/12*M*((2L^2) + c^2)) 
Kf = 0.0275;
Mh= 0.08;   %0.4
g = 9.81;
L1 = 0.3;
L2 = 0.1;
c = 0.0076;


L1 = 0.3;
L2 = 0.1;
M1 = 0.04;
M2 = 0.04;
G = 9.81;


A = [  0  1  0  0;
      -Mh*g*L1/Jteta   -c/Jteta  0  0;
       0  0  0  1;
       0  0  0  -c/Jphi];

B = [ 0  0  ;
      Kf*L1/Jteta 0;
      0  0 ;
      0 Kf*L2/Jphi];

C1 = [1 0 0 0];
C2 = [0 0 1 0];
C = [1 0 0 0;
     0 0 1 0];
 
D = C;

S1 = ss(A,B,C1,0);
S2 = ss(A,B,C2,0);

H1 = tf(S1);
H2 = tf(S2);
      
H11 = H1(1,1);      
H22 = H2(1,2);

% LQR

Q = eye(4);
Q(2,2) = 100;
Q(4,4) = 100;
R = eye(2);
R(1,1) = 0.1;

K = lqr(A,B,Q,R);

%% LQR integer pole
Az = [A, zeros(4,2)
      D, zeros(2,2)];
  
Bz = [B;
      zeros(2,2)];
  
Qz = [zeros(4,4), zeros(4,2)
      zeros(2,4), eye(2)];

Qz(5,5) = 1;   
Qz(6,6) = 1;
Rz = eye(2);

Kpi = lqr(Az,Bz,Qz,Rz);
K1 = Kpi(:,1:4);
K2 = Kpi(:,5:6);
      


