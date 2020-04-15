%SDRE Aero PEndulo

A = [  0  1  0  0;
      (Mw*g*Lw - Mh*g*L1)sinc(teta/pi)*pi/Jteta   -c/Jteta  0  0;
       0  0  0  1;
       0  0  0  -c/Jphi];

B = [ 0  0  ;
      Kf*L1*cos(phi)/Jteta 0;
      0  0 ;
      0 Kf*L2/Jphi];
  
  
  K = lqr(A,B