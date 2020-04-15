%%Teste model_Leis de Newton
%% Aero Pêndulo dual

Jteta = 0.005425;   % (1/12*M*((b^2) + c^2)) + M*L1^2 
Jphi = 0.00201;     % (1/12*M*((2L^2) + c^2)) 
Kf = 0.06; %0.037; 0.6     
% 0.0275;
Mh= 0.08;
Mw = 0.08;
Lw = 0.13;%0.13
g = 9.81;
L1 = 0.227;         
L2 = 0.08;  
c = 0.0076; %0.016
teta30 = 30*pi/180;
teta90 = 90*pi/180;

tetax = [10 30 50 70 90 110 130 150 170]*pi/180;
phix = [-80 -60 -40 -20 0 20 40 60 80]*pi/180;
L1 = 0.227;
L2 = 0.08;
M1 = 0.2;
M2 = 0.2;
G = 9.81;

% A = [  0  1  0  0;
%       (Mw*g*Lw - Mh*g*L1)*cos(teta)/(Jteta)   -c/Jteta  0  0;
%        0  0  0  1;
%        0  0  0  -c/Jphi];
 

B = [ 0  0  ;
      Kf*L1/Jteta 0;
      0  0 ;
      0 Kf*L2/Jphi];

C = [1 0 0 0;
     0 0 1 0];

D = C;

Bz = [B;
      zeros(2,2)];
  
Qz = [zeros(4,4), zeros(4,2)
      zeros(2,4), eye(2)];

Qz(5,5) = 20;   %20
Qz(6,6) = 5;   %10
Rz = eye(2);

Kpx = zeros(2*9,6*9);
Kp2 = zeros(2*9,6*9);

for i = 1:9
    Ax = [ 0 1 0 0;
       (Mw*g*Lw - Mh*g*L1)*cos(tetax(i))/(Jteta)  -c/Jteta  0  0;
       0  0  0  1;
       0  0  0  -c/Jphi];
   
    Azx = [Ax, zeros(4,2)
         D, zeros(2,2)];
   
   for j = 1:9
       
    Bx = [ 0  0  ;
      Kf*L1*cos(phix(j))/Jteta 0;
      0  0 ;
      0 Kf*L2/Jphi];
  
    Bzx = [Bx;
      zeros(2,2)];
   
    if(isempty(K1))
        pp00 =  1.792; pp01 = 0.05293; pp02 =  -4.347e-05;
        pp20 = -0.003316; pp21 = -2.997e-05; pp22 = 4.784e-07;
        pp40 =5.41e-07;
        p00 = 1.275; p01 = 0.001739; p02 = 2.001e-05;  
        p20 = -0.001137; p21 = 3.313e-06; p22 = 5.662e-08;  
        p40 = 1.397e-07;

        K1 = [0.345552654675701 0.207922866674390 1.52365480300072e-16 1.41750121698813e-17 4.4712 0;
              2.34353074450996e-17 1.34831465617473e-17 3.19016248053525 0.692343603732364 0 2.2361];
    end
    x = phix(j);
    y = tetax(i);
    K1(1,1) = pp00 + pp01*y + pp20*x^2 + pp02*y^2 + pp21*x^2*y + pp40*x^4 + pp22*x^2*y^2;
    K1(1,2) = p00 + p01*y + p20*x^2 + p02*y^2 + p21*x^2*y + p40*x^4 + p22*x^2*y^2;
  
    Kpx(i*2-1:i*2,6*j-5:6*j) = lqr(Azx,Bzx,Qz,Rz);
    Kp2(i*2-1:i*2,6*j-5:6*j) = K1;
   end
end

    
AL30 = [ 0 1 0 0;
       (Mw*g*Lw - Mh*g*L1)*cos(teta30)/(Jteta)  -c/Jteta  0  0;
       0  0  0  1;
       0  0  0  -c/Jphi];
   
AL30 = [ 0 1 0 0;
       (Mw*g*Lw - Mh*g*L1)*cos(teta90)/(Jteta)  -c/Jteta  0  0;
       0  0  0  1;
       0  0  0  -c/Jphi];
   


%% LQR integer pole
Az30 = [AL30, zeros(4,2)
      D, zeros(2,2)];

Az90 = [AL30, zeros(4,2)
      D, zeros(2,2)];

[Kpi30,S30] = lqr(Az30,Bz,Qz,Rz);
[Kpi90,S30] = lqr(Az90,Bz,Qz,Rz);

K1 =  Kpi90(:,1:4);
K2 =  Kpi90(:,5:6);

K1_30 = Kpi30(:,1:4);
K2_30 = Kpi30(:,5:6);
g1 = zeros(9,9);
g2 = zeros(9,9);
h1 = zeros(9,9);
h2 = zeros(9,9);
e1 = zeros(9,9);
e2 = zeros(9,9);
for n = 1:9
    for m=1:9
       g1(m,n) = Kpx(m*2-1,6*n-5);
       g2(m,n) = Kpx(m*2-1,6*n-4);
       h1(m,n) = Kp2(m*2-1,6*n-5);
       h2(m,n) = Kp2(m*2-1,6*n-4);
       e1(m,n) = (g1(m,n)-h1(m,n))^2;
       e2(m,n) = (g2(m,n)-h2(m,n))^2;
    end
end
EMQ1 = sqrt(sum(sum(e1)))/81;  
EMQ2 = sqrt(sum(sum(e2)))/81;  

teta3d = zeros(1,81);
phi3d = zeros(1,81);
g13d = zeros(1,81);
g23d = zeros(1,81);
for i1 =1:9
    teta3d(9*i1-8:9*i1) = tetax(i1)*180/pi;
    phi3d(9*i1-8:9*i1)  = phix*180/pi;
    g13d(9*i1-8:9*i1) = g1(i1,:);
    g23d(9*i1-8:9*i1) = g2(i1,:);
end


