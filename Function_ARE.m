function y=areklmn_reactor(teta,phi)
% STATE SPACE COEFFICIENT MATRICES
Jteta = 0.005425;    %(1/12*M*((b^2) + c^2)) + M*L1^2 
Jphi = 0.00201;   %(1/12*M*((2L^2) + c^2)) 
Kf = 0.056;      % 0.0275;
Mh= 0.08;
Mw = 0.039;
Lw = 0.15;
g = 9.81;
L1 = 0.227        
L2 = 0.08;
c = 0.0076;
L1 = 0.227;
L2 = 0.08;

A = [  0  1  0  0;
      (Mw*g*Lw - Mh*g*L1)*sin(teta)/Jteta   -c/Jteta  0  0;
       0  0  0  1;
       0  0  0  -c/Jphi];
B = [ 0  0  ;
      Kf*L1*cos(phi)/Jteta 0;
      0  0 ;
      0 Kf*L2/Jphi];
D = [1 0 0 0;
     0 0 1 0];
C = D;
% State weighting matrix and Control
Q=eye(4);
Q(1,1) = 20;
Q(3,3) = 10;
R=eye(2); % This value will be inverse of
A = [A, zeros(4,2)
          D, zeros(2,2)];
B = [B;
     zeros(2,2)];

Qz = [zeros(4,4), zeros(4,2)
      zeros(2,4), eye(2)];

Qz(5,5) = 20;   
Qz(6,6) = 1;

x10=rand([6,6]);

x0=[x10];
options = optimset('Display','iter','MaxFunEvals',50000000,'MaxIter' ,100,'TolFun',1e-6,'FunValCheck','on');
% fsolve to solve Algebric Riccati
figure(1)
x=fsolve(@areklmn_reactor,x0(:),options);

% Extracting elements from x to get
k=x(:);

% % Converting vectors into matrices
P=reshape(k,[6,6])

% Calculating Proportional term
K=inv(R)*B'*P
% Calculating Integral term
K1 = K(:,1:4);
K2 = K(:,5:6);

end

function y=areklmn_reactor(x,teta,phi)
Jteta = 0.005425;    %1.0348 ;  (1/12*M*((b^2) + c^2)) + M*L1^2 
Jphi = 0.00201;  %0.0451;  (1/12*M*((2L^2) + c^2)) 
Kf = 0.056;      
% 0.0275;
Mh= 0.08;
Mw = 0.039;
Lw = 0.15;
g = 9.81;
L1 = 0.227          %0.227;
L2 = 0.08;
c = 0.0076;
L1 = 0.227;
L2 = 0.08;

A = [  0  1  0  0;
      (Mw*g*Lw - Mh*g*L1)*sin(teta)/Jteta   -c/Jteta  0  0;
       0  0  0  1;
       0  0  0  -c/Jphi];
B = [ 0  0  ;
      Kf*L1*cos(phi)/Jteta 0;
      0  0 ;
      0 Kf*L2/Jphi];
  D = [1 0 0 0;
         0 0 1 0];
  
% State weighting matrix and Control
Q=eye(4);
Q(1,1) = 20;
Q(3,3) = 10;
R=eye(2); % This value will be inverse of
A = [A, zeros(4,2)
          D, zeros(2,2)];

    B = [B;
          zeros(2,2)];

    Qz = [zeros(4,4), zeros(4,2)
          zeros(2,4), eye(2)];

    Qz(5,5) = 20;   
    Qz(6,6) = 1;

% Extracting elements from x to get
P= reshape(x,[6,6])
% Converting vectors into matrices

% Four simultaneous nonlinear matrix
y= -A'*P - P*A - Qz + P*B*inv(R)*B'*P;
end