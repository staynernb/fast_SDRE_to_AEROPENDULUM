      REPS = 10000; minTime = Inf; nsum = 10;
      x = 20;
      y = 130;

      tic;
      for i=1:REPS
        tstart1 = tic;
        t1 = pp00 + pp01*y + pp20*x^2 + pp02*y^2 + pp21*x^2*y + pp40*x^4 + pp22*x^2*y^2;
        t2 = p00 + p01*y + p20*x^2 + p02*y^2 + p21*x^2*y + p40*x^4 + p22*x^2*y^2;
        telapsed1 = toc(tstart1);
        minTime1 = min(telapsed1,minTime);
      end
      averageTime1 = toc/REPS;
      
      tic;
      for i=1:REPS
        tstart2 = tic;
        Kpi30 = lqr(Az30,Bz,Qz,Rz);
        telapsed2 = toc(tstart2);
        minTime2 = min(telapsed2,minTime);
      end
      averageTime2 = toc/REPS;
      ganhoA = averageTime2/averageTime1;
      ganhoM = minTime2/minTime1;
      