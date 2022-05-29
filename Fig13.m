%% Fig.13A Perturbing u0 at the beginning of the closing phase: 
M=lyttle_model;
M.tmax=15; M.solve; 
M_period = M.t_open_to_close(end)-M.t_open_to_close(end-1);

% start from closing:
M.xinit = M.y_open_to_close(end,1:8); % start at close

M.tmax = 2*M_period; % run for two periods
M.solve

M3=lyttle_model;
eps=0.1;
M3.xinit = M.xinit;
M3.xinit(4)=M.xinit(4)+eps; M3.tmax=15; M3.solve;

M3_period = M3.t_open_to_close(end)-M3.t_open_to_close(end-1);

M3.tmax = 2*M3_period; % run for two periods
M3.solve

figure
set(gcf,'Position',[500 800 1200 400])
% Plot trajectories together, aligned from closing
%  plot brain variables
subplot(1,2,1)
plot(M.t,M.yext(:,1),'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on
plot(M.t,M.yext(:,2),'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(M.t,M.yext(:,3), 'Color',[0.9290 0.6940 0.1250],'LineWidth',2)

plot(M3.t,M3.yext(:,1),'--','Color',[0 0.4470 0.7410],'LineWidth',2)
plot(M3.t,M3.yext(:,2),'--','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(M3.t,M3.yext(:,3),'Color',[0.9290 0.6940 0.1250],'LineStyle','--','LineWidth',2)

plot([M3.t_close_to_open(1) M3.t_close_to_open(1)],[-0 1],'m','linewidth',2)
plot([M3.t_close_to_open(2) M3.t_close_to_open(2)],[-0 1],'m','linewidth',2)
title('Perturbing $u_0$ at the beginning of the closing phase leads to a phase delay','interpreter','latex','fontsize',24)
legend('$a_0$','$a_1$','$a_2$','Location','southwest','interpreter','latex')
set(gca,'FontSize',20)
xlim([M.t(1) M3.t(end)])
ylim([0, 1])
M.draw_wall_closing
xlabel('time','interpreter','latex','fontsize',20)

subplot(1,2,2)
plot(M.t,M.yext(:,4),'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on
plot(M.t,M.yext(:,5),'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(M.t,M.yext(:,6),'Color',[0.9290 0.6940 0.1250],'LineWidth',2)

plot(M3.t,M3.yext(:,4),'--','Color',[0 0.4470 0.7410],'LineWidth',2)
plot(M3.t,M3.yext(:,5),'--','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(M3.t,M3.yext(:,6),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',2)

plot([M3.t_close_to_open(1) M3.t_close_to_open(1)],[-0.4 1],'m','linewidth',2)
plot([M3.t_close_to_open(2) M3.t_close_to_open(2)],[-0.4 1],'m','linewidth',2)

xlabel('time','interpreter','latex','fontsize',20)
legend('$u_0$','$u_1$','$x_r$','Location','southwest','interpreter','latex')
set(gca,'FontSize',20)
xlim([M.t(1) M3.t(end)])
ylim([-0.4, 1])
M.draw_wall_closing


%% Fig.13B Perturbing u0 during the late closing phase 
M=lyttle_model;
M.tmax=15; M.solve; 
M_period = M.t_open_to_close(end)-M.t_open_to_close(end-1);

% unperturbed solution, start from closing
M.xinit = M.y_open_to_close(end,1:8); 
M.tmax = 2*M_period; % run for two periods
M.solve

% Perturb u0 at time 2.106
M3=lyttle_model;
M3.tinit=2.106;
ind=(abs(M.t-2.106)<0.0005);
M3.xinit = M.yext(ind,1:8);

eps=0.1;
M3.xinit(4)=M3.xinit(4)+eps; M3.tmax=M.tmax-M3.t0; M3.solve;


figure
set(gcf,'Position',[500 800 1200 400])
% Plot trajectories together, aligned from closing
%  plot brain variables
subplot(1,2,1)
plot(M.t,M.yext(:,1),'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on
plot(M.t,M.yext(:,2),'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(M.t,M.yext(:,3), 'Color',[0.9290 0.6940 0.1250],'LineWidth',2)

plot(M3.t,M3.yext(:,1),'--','Color',[0 0.4470 0.7410],'LineWidth',2)
plot(M3.t,M3.yext(:,2),'--','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(M3.t,M3.yext(:,3),'Color',[0.9290 0.6940 0.1250],'LineStyle','--','LineWidth',2)

plot([M3.t_close_to_open(1) M3.t_close_to_open(1)],[-0 1],'m','linewidth',2)
plot([M3.t_close_to_open(2) M3.t_close_to_open(2)],[-0 1],'m','linewidth',2)
title('Perturbing $u_0$ during the late closing phase leads to a phase advance','interpreter','latex','fontsize',24)
legend('$a_0$','$a_1$','$a_2$','Location','southwest','interpreter','latex')
set(gca,'FontSize',20)
xlim([M.t(1) M.tmax])
ylim([0, 1])
M.draw_wall_closing
xlabel('time','interpreter','latex','fontsize',20)

subplot(1,2,2)
plot(M.t,M.yext(:,4),'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on
plot(M.t,M.yext(:,5),'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(M.t,M.yext(:,6),'Color',[0.9290 0.6940 0.1250],'LineWidth',2)

plot(M3.t,M3.yext(:,4),'--','Color',[0 0.4470 0.7410],'LineWidth',2)
plot(M3.t,M3.yext(:,5),'--','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(M3.t,M3.yext(:,6),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',2)

plot([M3.t_close_to_open(1) M3.t_close_to_open(1)],[-0.4 1],'m','linewidth',2)
plot([M3.t_close_to_open(2) M3.t_close_to_open(2)],[-0.4 1],'m','linewidth',2)

xlabel('time','interpreter','latex','fontsize',20)
legend('$u_0$','$u_1$','$x_r$','Location','southwest','interpreter','latex')
set(gca,'FontSize',20)
xlim([M.t(1) M.tmax])
ylim([-0.4, 1])
M.draw_wall_closing

%% Consider perturbing  u1 at the beginning of the closing phase Fig 13C
M=lyttle_model;
M.tmax=15; M.solve; 
M_period = M.t_open_to_close(end)-M.t_open_to_close(end-1);


% start from closing:
M.xinit = M.y_open_to_close(end,1:8); % start at close

M.tmax = 2*M_period; % run for two periods
M.solve


M3=lyttle_model;
eps=0.1;
M3.xinit = M.xinit;
M3.xinit(5)=M.xinit(5)+eps; M3.tmax=15; M3.solve;

M3_period = M3.t_open_to_close(end)-M3.t_open_to_close(end-1);

M3.tmax = 2*M3_period; % run for two periods
M3.solve

figure
set(gcf,'Position',[500 800 1200 400])
% Plot trajectories together, aligned from closing
%  plot brain variables
subplot(1,2,1)
plot(M.t,M.yext(:,1),'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on
plot(M.t,M.yext(:,2),'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(M.t,M.yext(:,3), 'Color',[0.9290 0.6940 0.1250],'LineWidth',2)

plot(M3.t,M3.yext(:,1),'--','Color',[0 0.4470 0.7410],'LineWidth',2)
plot(M3.t,M3.yext(:,2),'--','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(M3.t,M3.yext(:,3),'Color',[0.9290 0.6940 0.1250],'LineStyle','--','LineWidth',2)

plot([M3.t_close_to_open(1) M3.t_close_to_open(1)],[-0 1],'m','linewidth',2)
plot([M3.t_close_to_open(2) M3.t_close_to_open(2)],[-0 1],'m','linewidth',2)

legend('$a_0$','$a_1$','$a_2$','Location','southwest','interpreter','latex')
set(gca,'FontSize',20)
xlim([M.t(1) M3.t(end)])
ylim([0, 1])
M.draw_wall_closing
xlabel('time','interpreter','latex','fontsize',20)
title('Perturbing $u_1$ at the beginning of the closing phase leads to a phase advance','interpreter','latex','fontsize',24)

subplot(1,2,2)
plot(M.t,M.yext(:,4),'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on
plot(M.t,M.yext(:,5),'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(M.t,M.yext(:,6),'Color',[0.9290 0.6940 0.1250],'LineWidth',2)

plot(M3.t,M3.yext(:,4),'--','Color',[0 0.4470 0.7410],'LineWidth',2)
plot(M3.t,M3.yext(:,5),'--','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(M3.t,M3.yext(:,6),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',2)

plot([M3.t_close_to_open(1) M3.t_close_to_open(1)],[-0.4 1],'m','linewidth',2)
plot([M3.t_close_to_open(2) M3.t_close_to_open(2)],[-0.4 1],'m','linewidth',2)

xlabel('time','interpreter','latex','fontsize',20)
legend('$u_0$','$u_1$','$x_r$','Location','southwest','interpreter','latex')
set(gca,'FontSize',20)
xlim([M.t(1) M3.t(end)])
ylim([-0.4, 1])
M.draw_wall_closing

%% Perturbing u1 during the late closing phase Fig 13D
M=lyttle_model;
M.tmax=15; M.solve; 
M_period = M.t_open_to_close(end)-M.t_open_to_close(end-1);

% unperturbed solution, start from closing
M.xinit = M.y_open_to_close(end,1:8); 
M.tmax = 2*M_period; % run for two periods
M.solve

% Perturb u0 at time 2.106
M3=lyttle_model;
M3.tinit=2.106;
ind=(abs(M.t-2.106)<0.0005);
M3.xinit = M.yext(ind,1:8);

eps=0.1;
M3.xinit(5)=M3.xinit(5)+eps; M3.tmax=M.tmax-M3.t0; M3.solve;


figure
set(gcf,'Position',[500 800 1200 400])
% Plot trajectories together, aligned from closing
%  plot brain variables
subplot(1,2,1)
plot(M.t,M.yext(:,1),'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on
plot(M.t,M.yext(:,2),'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(M.t,M.yext(:,3), 'Color',[0.9290 0.6940 0.1250],'LineWidth',2)

plot(M3.t,M3.yext(:,1),'--','Color',[0 0.4470 0.7410],'LineWidth',2)
plot(M3.t,M3.yext(:,2),'--','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(M3.t,M3.yext(:,3),'Color',[0.9290 0.6940 0.1250],'LineStyle','--','LineWidth',2)

plot([M3.t_close_to_open(1) M3.t_close_to_open(1)],[-0 1],'m','linewidth',2)
plot([M3.t_close_to_open(2) M3.t_close_to_open(2)],[-0 1],'m','linewidth',2)
title('Perturbing $u_1$ during the late closing phase leads to a phase delay','interpreter','latex','fontsize',24)
legend('$a_0$','$a_1$','$a_2$','Location','southwest','interpreter','latex')
set(gca,'FontSize',20)
xlim([M.t(1) M.tmax])
ylim([0, 1])
M.draw_wall_closing
xlabel('time','interpreter','latex','fontsize',20)

subplot(1,2,2)
plot(M.t,M.yext(:,4),'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on
plot(M.t,M.yext(:,5),'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(M.t,M.yext(:,6),'Color',[0.9290 0.6940 0.1250],'LineWidth',2)

plot(M3.t,M3.yext(:,4),'--','Color',[0 0.4470 0.7410],'LineWidth',2)
plot(M3.t,M3.yext(:,5),'--','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(M3.t,M3.yext(:,6),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',2)

plot([M3.t_close_to_open(1) M3.t_close_to_open(1)],[-0.4 1],'m','linewidth',2)
plot([M3.t_close_to_open(2) M3.t_close_to_open(2)],[-0.4 1],'m','linewidth',2)

xlabel('time','interpreter','latex','fontsize',20)
legend('$u_0$','$u_1$','$x_r$','Location','southwest','interpreter','latex')
set(gca,'FontSize',20)
xlim([M.t(1) M.tmax])
ylim([-0.4, 1])
M.draw_wall_closing
