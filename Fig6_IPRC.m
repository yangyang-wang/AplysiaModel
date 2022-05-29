% Generate and plot IPRC (Fig.6) and compute the global linear change in period T1 by running lyttle_iPRC 
eps1 = 1e-4;
fsw=0.01;
pl = true; % plot

[~,~,T1] = lyttle_iPRC(eps1,fsw,pl);