% This file computes and plots iPRC. 
function [t_iprc,x_iprc,T1] = lyttle_iPRC(eps1,fsw,pl)

if nargin<1
    eps1 = 1e-4;
end
if nargin<2
    fsw=0.01;
end
if nargin<3
    pl = false;
end


% For different eps values, find y_open_to_close, y_close_to_open and the total closing time. 

if eps1 == 0.0001 && fsw == 0.01
    model0 = lyttle_model;
    model0.solve;
    T0=model0.tmax;
    
    model=lyttle_model('xinit', model0.y_open_to_close(end,1:8)); % start at the closing 
    model.solve
else
% Need to run model0 for the transients to decay non-default eps1 and fsw values 

model0 = lyttle_model('eps1',eps1,'tmax',20); % new eps1
model0.xinit(8)=fsw; % update fsw 
model0.solve;
% restart at y_open_to_close
model1 = lyttle_model('eps1',eps1,'tmax',10,'xinit',model0.y_open_to_close(end,1:8)); % start at y_open_to_close
model1.solve
T0 = model1.findPeriod(4.5,6.5);

model = lyttle_model('eps1',eps1,'tmax',T0,'xinit',model1.y_open_to_close(end,1:8));
model.solve
end

% T0=model.tmax; % Period of LC solution of the Aplysia model when initial value is xinit. 

% first use find_prc_monodromy to find the initial cond for prc
M = find_prc_monodromy(model);
[V, D]=eig(M');
ind_e1=(abs(diag(D))==max(abs(diag(D)))); %find index where the eigenvalue is 1
zinit=V(:,ind_e1);
dummy=0;
f10=model.lyttle_ODE(dummy, model.yext(end,:),model.checkdomain(model.yext(end,:)));
f10=f10(1:6);
rescale=f10'*zinit;
prcinit=zinit'/(rescale);

model.find_prc(prcinit); % Use find_prc function to compute iPRC
% model.plot_prc;          % plot the results

t_iprc=model.prct;
x_iprc=model.prc;

% Below is to estimate relative changes in time in response to static perturbation
% applied to the seaweed: T1=int_{closing time} Z(t)\part F/\part eps (gamma(t),0)dt = int_{closing time}z_6(t)*(1)

a1=model.yext(:,2);
a2=model.yext(:,3);
ind_close=(a1+a2>0.5);
int=-trapz(model.prct(flipud(ind_close)),model.prc(flipud(ind_close),6)/model.br);% time is backward,so the integral has the opposite sign 
T1=-int; 
% nu=T1/T0;
% disp('T1 is')
% disp(T1)
% disp('nu is')
% disp(nu)

%% Plot the results
if pl
    
% figure
% set(gcf,'Position',[50 800 800 800])
% subplot(2,2,1)
% plot(model.t,model.yext(:,1:3),'linewidth',2)
% xlim([0 T0])
% ylim([-0.01 1.1])
% set(gca,'FontSize',18)
% legend('$a_0$','$a_1$','$a_2$','Location','southwest','interpreter','latex')
% model.draw_wall_closing           
% 
% subplot(2,2,3)
% plot(model.prct,model.prc(:,1:3),'linewidth',2)
% xlim([0 T0])
% ylim([-0.01 5.3e+4])
% set(gca,'FontSize',18)
% legend('$z_{a_0}$','$z_{a_1}$','$z_{a_2}$','interpreter','latex','fontsize',25)
% model.draw_wall_closing
% 
% subplot(2,2,2)
% plot(model.t,model.yext(:,4:6),'linewidth',2)
% xlim([0 T0])
% ylim([-0.01 1.1])
% set(gca,'FontSize',18)
% legend('$u_0$','$u_1$','$x_r$','Location','southwest','interpreter','latex')
% model.draw_wall_closing         
% 
% subplot(2,2,4)
% plot(model.prct,model.prc(:,4:6),'linewidth',2)
% ylim([-5 5])
% xlim([0 T0])
% xlabel('time','interpreter','latex','fontsize',25)
% legend('$z_{u_0}$','$z_{u_1}$','$z_{x_r}$','interpreter','latex','fontsize',25)
% set(gca,'FontSize',18)
% model.draw_wall_closing

figure
subplot(2,1,1)
plot(model.prct,model.prc(:,1:3),'linewidth',2)
xlim([0 T0])
ylim([-0.01 5.3e+4])
title('IPRC','Interpreter','latex','FontWeight','normal','Fontsize',20)
set(gca,'FontSize',18)
legend('$z_{a_0}$','$z_{a_1}$','$z_{a_2}$','interpreter','latex','fontsize',25)
model.draw_wall_closing

subplot(2,1,2)
plot(model.prct,model.prc(:,4:6),'linewidth',2)
hold on; plot([0 T0],[0 0],'k:','linewidth',2)
ylim([-3 6])
xlim([0 T0])
xlabel('time','interpreter','latex','fontsize',25)
legend('$z_{u_0}$','$z_{u_1}$','$z_{x_r}$','interpreter','latex','fontsize',25)
set(gca,'FontSize',18)
model.draw_wall_closing


end
end