tmax=200;
alpha=.2;
beta=2;
yinit=[alpha,-1];
[t,y]=LC_one_hard_bdy_one_switch(yinit,tmax,alpha,beta);
figure
plot(y(:,1),y(:,2))
grid on
shg
if ~prod(~(abs(axis)>100))
    xlim([-10 10])
end