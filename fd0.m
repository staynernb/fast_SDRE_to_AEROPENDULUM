%ARE SOlution


function dt = are_s
tf = 25;
t0 = 0;
x0 = 0;
options = odeset('Abstol', 1e-6, 'Reltol',1e-6)
[t y] = ode45(@fdo, [t0 tf], x0, options);
plot(t,y);

function dP = fdo(P)
persistent A;
if isempty(A)
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
    A = [  0  1  0  0;
          (Mw*g*Lw - Mh*g*L1)/Jteta   -c/Jteta  0  0;
           0  0  0  1;
           0  0  0  -c/Jphi];

    B = [ 0  0  ;
          Kf*L1/Jteta 0;
          0  0 ;
          0 Kf*L2/Jphi];
    D = [1 0 0 0;
         0 0 1 0];

    A = [A, zeros(4,2)
          D, zeros(2,2)];

    B = [B;
          zeros(2,2)];

    Qz = [zeros(4,4), zeros(4,2)
          zeros(2,4), eye(2)];

    Qz(5,5) = 20;   
    Qz(6,6) = 1;
    R= eye(2);
end

dP = -A'*P - P*A - Qz + P*B*inv(R)*B'*P;