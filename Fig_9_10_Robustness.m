% Figs 9 and 10
% Compute Robustness and Task-fitness as sensory feedback gains or muscle strengths vary 


%% Fig.9 top
eps1_range = (-6:0.1:-2);
% eps1_range = (-4.8:0.1:-3); % zoom in


fsw_range = 0.01;
eps = 0.001; % fsw -> fsw+eps

T0=nan(length(eps1_range),length(fsw_range));
delta_x0=T0;
Teps=T0;
delta_xeps=T0;
Seps=T0;
S0=T0;
dx1_over_dx0 = T0;
T1_over_T0 = T0;

tic;

for i=1:length(eps1_range)
    for j=1:length(fsw_range)
        
        eps1=10^eps1_range(i);
        fsw = fsw_range(j);
        
        % no perturbation, with parameters eps1 and fsw
        model1 = lyttle_model('tmax',20,'eps1',eps1); % assign new eps1 value
        model1.xinit(8)=fsw; % update fsw
        model1.solve;
        model = lyttle_model('tmax',10,'eps1',eps1,'xinit',model1.y_open_to_close(end,1:8)); % start at y_open_to_close
        model.solve
%         T0(i,j)=model.findPeriod(4.5,6.5);
        disp(['i is ' num2str(i) ' and j is ' num2str(j)])
        
        if isempty(model.t_open_to_close)
            continue
        end
            
        T0(i,j)=model.t_open_to_close(1);
        
        if model.t_open_to_close(1)<0.1 && length(model.t_open_to_close)>1
            T0(i,j)=model.t_open_to_close(2); % in case model takes the initial starting time as the first t_open_to_close
        elseif model.t_open_to_close(1)<0.1 && length(model.t_open_to_close)==1
            T0(i,j) = nan; 
            continue
        end
        
        delta_x0(i,j) = model.y_close_to_open(1,6) - model.xinit(1,6);   % netchange of grasper position over the closing phase
        S0(i,j) = delta_x0(i,j)/T0(i,j);
        
        % with perturbation
        model1.xinit(8) = fsw+eps;
        model1.solve;
        model_pert = lyttle_model('tmax',10,'eps1',eps1,'xinit',model1.y_open_to_close(end,1:8)); % start at y_open_to_close
        model_pert.solve
%         Teps(i,j)=model_pert.findPeriod(4.5,6.5);
        Teps(i,j)=model_pert.t_open_to_close(1);
        delta_xeps(i,j) = model_pert.y_close_to_open(1,6) - model_pert.xinit(1,6);
        Seps(i,j) = delta_xeps(i,j)/Teps(i,j);
        
        % compute T1 and delta_x1
%         T1 = lyttle_iPRC(eps1,fsw); % T1=8.0778; when eps1=1e-4, fsw=0.01
        T1 = (Teps(i,j)-T0(i,j))/eps;
        delta_x1 = (delta_xeps(i,j) - delta_x0(i,j))/eps;
        
        dx1_over_dx0(i,j) = delta_x1/delta_x0(i,j);
        T1_over_T0(i,j) = T1/T0(i,j);                     
    end
end

robustness = fsw*(Seps-S0)./(eps*S0);  % Equation 2.5

toc 

% save intake-rate.mat eps1_range fsw_range task_fitness_1storder task_fitness T0 Teps S0 Seps delta_x0 delta_xeps

% Robustness/sensitivity to pertbation on fsw
figure
yyaxis left
plot(eps1_range,robustness,'k','linewidth',2);
xlabel('$\rm log10(\varepsilon_0)$','interpreter','latex','fontsize',30)
ylabel('$\frac{d S_\varepsilon}{d\varepsilon}/\varepsilon_0$','interpreter','latex','fontsize',30,'rot',0)
title('robustness' ,'FontWeight','normal','interpreter','latex','fontsize',25)
% xlim([eps1_range(1)-0.001 eps1_range(end)])

yyaxis right
plot(eps1_range,-S0,'b','linewidth',2) % -S0 = -deltaX/T 
% xlim([eps1_range(1)-0.001 eps1_range(end)])
ylabel('$S_0$','Interpreter', 'Latex','Fontsize',15,'rot',0);
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
set(ax,'FontSize',18)

%% Fig.9 middle
eps2_range = (-6:0.1:-2);
fsw_range = 0.01;
eps = 0.001; % fsw -> fsw+eps

T0=nan(length(eps2_range),length(fsw_range));
delta_x0=T0;
Teps=T0;
delta_xeps=T0;
Seps=T0;
S0=T0;
dx1_over_dx0 = T0;
T1_over_T0 = T0;

tic;

for i=1:length(eps2_range)
    for j=1:length(fsw_range)
        
        eps2=10^eps2_range(i);
        fsw = fsw_range(j);
        
        % no perturbation, with parameters eps1 and fsw
        model1 = lyttle_model('tmax',20,'eps2',eps2); % assign new eps1 value
        model1.xinit(8)=fsw; % update fsw
        model1.solve;
        model = lyttle_model('tmax',10,'eps2',eps2,'xinit',model1.y_open_to_close(end,1:8)); % start at y_open_to_close
        model.solve
%         T0(i,j)=model.findPeriod(4.5,6.5);
        disp(['i is ' num2str(i) ' and j is ' num2str(j)])
        
        if isempty(model.t_open_to_close)
            continue
        end
            
        T0(i,j)=model.t_open_to_close(1);
        
        if model.t_open_to_close(1)<0.1 && length(model.t_open_to_close)>1
            T0(i,j)=model.t_open_to_close(2); % in case model takes the initial starting time as the first t_open_to_close
        elseif model.t_open_to_close(1)<0.1 && length(model.t_open_to_close)==1
            T0(i,j) = nan; 
            continue
        end
        
        delta_x0(i,j) = model.y_close_to_open(1,6) - model.xinit(1,6);   % netchange of grasper position over the closing phase
        S0(i,j) = delta_x0(i,j)/T0(i,j);
        
        % with perturbation
        model1.xinit(8) = fsw+eps;
        model1.solve;
        model_pert = lyttle_model('tmax',10,'eps2',eps2,'xinit',model1.y_open_to_close(end,1:8)); % start at y_open_to_close
        model_pert.solve
%         Teps(i,j)=model_pert.findPeriod(4.5,6.5);
        Teps(i,j)=model_pert.t_open_to_close(1);
        delta_xeps(i,j) = model_pert.y_close_to_open(1,6) - model_pert.xinit(1,6);
        Seps(i,j) = delta_xeps(i,j)/Teps(i,j);
        
        % compute T1 and delta_x1
%         T1 = lyttle_iPRC(eps1,fsw); % T1=8.0778; when eps1=1e-4, fsw=0.01
        T1 = (Teps(i,j)-T0(i,j))/eps;
        delta_x1 = (delta_xeps(i,j) - delta_x0(i,j))/eps;
        
        dx1_over_dx0(i,j) = delta_x1/delta_x0(i,j);
        T1_over_T0(i,j) = T1/T0(i,j);                     
    end
end

robustness = fsw*(Seps-S0)./(eps*S0);  

toc 

% Robustness/sensitivity to pertbation on fsw
figure
yyaxis left
plot(eps2_range,robustness,'k','linewidth',2);
xlabel('$\rm log10(\varepsilon_1)$','interpreter','latex','fontsize',30)
ylabel('$\frac{d S_\varepsilon}{d\varepsilon}/S_0$','interpreter','latex','fontsize',30,'rot',0)
title('robustness' ,'FontWeight','normal','interpreter','latex','fontsize',25)
% xlim([eps1_range(1)-0.001 eps1_range(end)])

yyaxis right
plot(eps2_range,-S0,'b','linewidth',2)
% xlim([eps1_range(1)-0.001 eps1_range(end)])
ylabel('$S_0$','Interpreter', 'Latex','Fontsize',15,'rot',0);
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
set(ax,'FontSize',18)

%% Fig.9 bottom
% eps3_range = linspace(5e-5,1e-2,20);
eps3_range = (-6:0.1:-2);
fsw_range = 0.01;
eps = 0.001; % fsw -> fsw+eps

T0=nan(length(eps3_range),length(fsw_range));
delta_x0=T0;
Teps=T0;
delta_xeps=T0;
Seps=T0;
S0=T0;
dx1_over_dx0 = T0;
T1_over_T0 = T0;

tic;

for i=1:length(eps3_range)
    for j=1:length(fsw_range)
        
        eps3=10^eps3_range(i);
        fsw = fsw_range(j);
        
        % no perturbation, with parameters eps1 and fsw
        model1 = lyttle_model('tmax',20,'eps3',eps3); % assign new eps1 value
        model1.xinit(8)=fsw; % update fsw
        model1.solve;
        model = lyttle_model('tmax',10,'eps3',eps3,'xinit',model1.y_open_to_close(end,1:8)); % start at y_open_to_close
        model.solve
%         T0(i,j)=model.findPeriod(4.5,6.5);
        disp(['i is ' num2str(i) ' and j is ' num2str(j)])
        
        if isempty(model.t_open_to_close)
            continue
        end
            
        T0(i,j)=model.t_open_to_close(1);
        
        if model.t_open_to_close(1)<0.1 && length(model.t_open_to_close)>1
            T0(i,j)=model.t_open_to_close(2); % in case model takes the initial starting time as the first t_open_to_close
        elseif model.t_open_to_close(1)<0.1 && length(model.t_open_to_close)==1
            T0(i,j) = nan; 
            continue
        end
        
        delta_x0(i,j) = model.y_close_to_open(1,6) - model.xinit(1,6);   % netchange of grasper position over the closing phase
        S0(i,j) = delta_x0(i,j)/T0(i,j);
        
        % with perturbation
        model1.xinit(8) = fsw+eps;
        model1.solve;
        model_pert = lyttle_model('tmax',10,'eps3',eps3,'xinit',model1.y_open_to_close(end,1:8)); % start at y_open_to_close
        model_pert.solve
%         Teps(i,j)=model_pert.findPeriod(4.5,6.5);
        Teps(i,j)=model_pert.t_open_to_close(1);
        delta_xeps(i,j) = model_pert.y_close_to_open(1,6) - model_pert.xinit(1,6);
        Seps(i,j) = delta_xeps(i,j)/Teps(i,j);
        
        % compute T1 and delta_x1
%         T1 = lyttle_iPRC(eps1,fsw); % T1=8.0778; when eps1=1e-4, fsw=0.01
        T1 = (Teps(i,j)-T0(i,j))/eps;
        delta_x1 = (delta_xeps(i,j) - delta_x0(i,j))/eps;
        
        dx1_over_dx0(i,j) = delta_x1/delta_x0(i,j);
        T1_over_T0(i,j) = T1/T0(i,j);                     
    end
end

robustness = fsw*(Seps-S0)./(eps*S0);

toc 

% save intake-rate.mat eps1_range fsw_range task_fitness_1storder task_fitness T0 Teps S0 Seps delta_x0 delta_xeps

% Robustness/sensitivity to pertbation on fsw
figure
yyaxis left
plot(eps3_range,robustness,'k','linewidth',2);
xlabel('$\rm log10(\varepsilon_2)$','interpreter','latex','fontsize',30)
ylabel('$\frac{d S_\varepsilon}{d\varepsilon}$','interpreter','latex','fontsize',30,'rot',0)
title('robustness' ,'FontWeight','normal','interpreter','latex','fontsize',25)
% xlim([eps1_range(1)-0.001 eps1_range(end)])

yyaxis right
plot(eps3_range,-S0,'b','linewidth',2)
% xlim([eps1_range(1)-0.001 eps1_range(end)])
ylabel('$S_0$','Interpreter', 'Latex','Fontsize',15,'rot',0);
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
set(ax,'FontSize',18)


%% Fig.10 Left 
%explore sensitivity to mechanical changes: e.g., k0 = -1

k0_range = linspace(0.5,3,20);
fsw_range = 0.01;
eps = 0.001; % fsw -> fsw+eps

T0=nan(length(k0_range),length(fsw_range));
T1=T0;
delta_x1=T0;
delta_x0=T0;
Teps=T0;
delta_xeps=T0;
Seps=T0;
S0=T0;
dx1_over_dx0 = T0;
T1_over_T0 = T0;

T0_close =T0;

tic;

for i=1:length(k0_range)
    for j=1:length(fsw_range)
        
        k0=k0_range(i);
        fsw = fsw_range(j);
        
        % no perturbation, with parameters k0 and fsw
        model1 = lyttle_model('tmax',20); % assign new k0 value
        model1.k0 = k0;
        model1.xinit(8)=fsw; % update fsw
        model1.solve;
        model = lyttle_model('tmax',10,'xinit',model1.y_open_to_close(end,1:8)); % start at y_open_to_close
        model.k0 = k0;
        model.solve
%         T0(i,j)=model.findPeriod(4.5,6.5);
        disp(['i is ' num2str(i) ' and j is ' num2str(j)])
        
        if isempty(model.t_open_to_close)
            continue
        end
            
        T0(i,j)=model.t_open_to_close(1);
        T0_close(i,j)=model.t_close_to_open(1); 
        
        if model.t_open_to_close(1)<0.1 && length(model.t_open_to_close)>1
            T0(i,j)=model.t_open_to_close(2); % in case model takes the initial starting time as the first t_open_to_close
        elseif model.t_open_to_close(1)<0.1 && length(model.t_open_to_close)==1
            T0(i,j) = nan; 
            continue
        end
        
        delta_x0(i,j) = model.y_close_to_open(1,6) - model.xinit(1,6);   % netchange of grasper position over the closing phase
        S0(i,j) = delta_x0(i,j)/T0(i,j);
        
        % with perturbation
        model1.xinit(8) = fsw+eps;
        model1.solve;
        model_pert = lyttle_model('tmax',10,'xinit',model1.y_open_to_close(end,1:8)); % start at y_open_to_close
        model_pert.k0 = k0;
        model_pert.solve
%         Teps(i,j)=model_pert.findPeriod(4.5,6.5);
        Teps(i,j)=model_pert.t_open_to_close(1);
        delta_xeps(i,j) = model_pert.y_close_to_open(1,6) - model_pert.xinit(1,6);
        Seps(i,j) = delta_xeps(i,j)/Teps(i,j);
        
        % compute T1 and delta_x1
%         T1 = lyttle_iPRC(eps1,fsw); % T1=8.0778; when eps1=1e-4, fsw=0.01
        T1(i,j) = (Teps(i,j)-T0(i,j))/eps;
        delta_x1(i,j) = (delta_xeps(i,j) - delta_x0(i,j))/eps;
        
        dx1_over_dx0(i,j) = delta_x1(i,j)/delta_x0(i,j);
        T1_over_T0(i,j) = T1(i,j)/T0(i,j);                     
    end
end

robustness = fsw*(Seps-S0)./(eps*S0);  % (Seps-S0)/(eps*S0)
ind = (abs(robustness)>1);
robustness(ind)=nan;

toc 


% Robustness/sensitivity to pertbation on fsw

figure
subplot(4,1,1)
yyaxis left
plot(k0_range,robustness,'k','linewidth',2);
xlabel('$k_0$','interpreter','latex','fontsize',30)
ylabel('robustness','interpreter','latex','fontsize',18)
xlim([k0_range(1)-0.001 k0_range(end)])
ylim([-0.02 -0.006])

yyaxis right
plot(k0_range,-S0,'b','linewidth',2)
xlim([k0_range(1)-0.001 k0_range(end)])
ylabel('$S_0$','Interpreter', 'Latex','Fontsize',15,'rot',0);
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
set(ax,'FontSize',18)
ylim([0.05 0.12])

subplot(4,1,2)
plot(k0_range,T0,'linewidth',2)
hold on
plot(k0_range,T1,'linewidth',2)
% plot(k0_range,T1_over_T0,'linewidth',2)
% plot(k0_range,T0_close,'linewidth',2)
legend('$T_0$','$T_1$','interpreter','latex','fontsize',18)
xlabel('$k_0$','interpreter','latex','fontsize',30)
set(gca,'FontSize',20)
xlim([k0_range(1)-0.001 k0_range(end)])
ylim([0 13])

subplot(4,1,3)
plot(k0_range,-delta_x0,'linewidth',2)
hold on
plot(k0_range,-delta_x1,'linewidth',2)
% plot(k0_range,dx1_over_dx0,'linewidth',2)
legend('$-\Delta x_{r,0}$','$-\Delta x_{r,1}$','interpreter','latex','fontsize',18)
xlabel('$k_0$','interpreter','latex','fontsize',30)
set(gca,'FontSize',20)
xlim([k0_range(1)-0.001 k0_range(end)])
ylim([0 1])

subplot(4,1,4)
plot(k0_range,T1_over_T0,'linewidth',2)
hold on
plot(k0_range,dx1_over_dx0,'linewidth',2)
legend('$T_1/T_0$','$\Delta x_{r,1}/\Delta x_{r,0}$','interpreter','latex','fontsize',18)
xlabel('$k_0$','interpreter','latex','fontsize',30)
set(gca,'FontSize',20)
xlim([k0_range(1)-0.001 k0_range(end)])
ylim([0 2.9])


%% Fig.10 Right
%explore sensitivity to mechanical changes: e.g., k1 = -1

k1_range = linspace(-3,-0.5,20);
fsw_range = 0.01;
eps = 0.001; % fsw -> fsw+eps

T0=nan(length(k1_range),length(fsw_range));
T1=T0;
delta_x1=T0;
delta_x0=T0;
Teps=T0;
delta_xeps=T0;
Seps=T0;
S0=T0;
dx1_over_dx0 = T0;
T1_over_T0 = T0;

T0_close = T0;

tic;

for i=1:length(k1_range)
    for j=1:length(fsw_range)
        
        k1=k1_range(i);
        fsw = fsw_range(j);
        
        % no perturbation, with parameters k1 and fsw
        model1 = lyttle_model('tmax',20); % assign new k1 value
        model1.k1 = k1;
        model1.xinit(8)=fsw; % update fsw
        model1.solve;
        model = lyttle_model('tmax',10,'xinit',model1.y_open_to_close(end,1:8)); % start at y_open_to_close
        model.k1 = k1;
        model.solve
%         T0(i,j)=model.findPeriod(4.5,6.5);
        disp(['i is ' num2str(i) ' and j is ' num2str(j)])
        
        if isempty(model.t_open_to_close)
            continue
        end
            
        T0(i,j)=model.t_open_to_close(1);
        
        T0_close(i,j)=model.t_close_to_open(1); 
        
        if model.t_open_to_close(1)<0.1 && length(model.t_open_to_close)>1
            T0(i,j)=model.t_open_to_close(2); % in case model takes the initial starting time as the first t_open_to_close
        elseif model.t_open_to_close(1)<0.1 && length(model.t_open_to_close)==1
            T0(i,j) = nan; 
            continue
        end
        
        delta_x0(i,j) = model.y_close_to_open(1,6) - model.xinit(1,6);   % netchange of grasper position over the closing phase
        S0(i,j) = delta_x0(i,j)/T0(i,j);
        
        % with perturbation
        model1.xinit(8) = fsw+eps;
        model1.solve;
        model_pert = lyttle_model('tmax',10,'xinit',model1.y_open_to_close(end,1:8)); % start at y_open_to_close
        model_pert.k1 = k1;
        model_pert.solve
%         Teps(i,j)=model_pert.findPeriod(4.5,6.5);
        Teps(i,j)=model_pert.t_open_to_close(1);
        delta_xeps(i,j) = model_pert.y_close_to_open(1,6) - model_pert.xinit(1,6);
        Seps(i,j) = delta_xeps(i,j)/Teps(i,j);
        
        % compute T1 and delta_x1
%         T1 = lyttle_iPRC(eps1,fsw); % T1=8.0778; when eps1=1e-4, fsw=0.01
        T1(i,j) = (Teps(i,j)-T0(i,j))/eps;
        delta_x1(i,j) = (delta_xeps(i,j) - delta_x0(i,j))/eps;
        
        dx1_over_dx0(i,j) = delta_x1(i,j)/delta_x0(i,j);
        T1_over_T0(i,j) = T1(i,j)/T0(i,j);                     
    end
end

robustness = fsw*(Seps-S0)./(eps*S0);  % (Seps-S0)/(eps*S0)
ind = (abs(robustness)>1);
robustness(ind)=nan;
toc 

% Robustness/sensitivity to pertbation on fsw
figure
subplot(4,1,1)
yyaxis left

plot(abs(k1_range),robustness,'k','linewidth',2);
xlabel('$|k_1|$','interpreter','latex','fontsize',30)
ylabel('Robustness','interpreter','latex','fontsize',18)
% xlim([eps1_range(1)-0.001 eps1_range(end)])
ylim([-0.02 -0.006])

yyaxis right
plot(abs(k1_range),-S0,'b','linewidth',2)
% xlim([eps1_range(1)-0.001 eps1_range(end)])
ylabel('$S_0$','Interpreter', 'Latex','Fontsize',15,'rot',0);
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
set(ax,'FontSize',18)
ylim([0.05 0.12])

subplot(4,1,2)
plot(abs(k1_range),T0,'linewidth',2)
hold on
plot(abs(k1_range),T1,'linewidth',2)
% plot(abs(k1_range),T1_over_T0,'linewidth',2)
% plot(abs(k1_range),T0_close,'linewidth',2)
% legend('$T_0$','$T_1$','$\frac{T1}{T0}$','$T_{0\rm clo}$','interpreter','latex','fontsize',18)
legend('$T_0$','$T_1$','interpreter','latex','fontsize',18)
xlabel('$|k_1|$','interpreter','latex','fontsize',30)
set(gca,'FontSize',20)
ylim([0 13])

subplot(4,1,3)
plot(abs(k1_range),-delta_x0,'linewidth',2)
hold on
plot(abs(k1_range),-delta_x1,'linewidth',2)
% plot(abs(k1_range),dx1_over_dx0,'linewidth',2)
legend('$-\Delta x_{r,0}$','$-\Delta x_{r,1}$','interpreter','latex','fontsize',15)
xlabel('$|k_1|$','interpreter','latex','fontsize',30)
set(gca,'FontSize',20)
ylim([0 1])

subplot(4,1,4)
plot(abs(k1_range),T1_over_T0,'linewidth',2)
hold on
plot(abs(k1_range),dx1_over_dx0,'linewidth',2)
legend('$T_1/T_0$','$\Delta x_{r,1}/\Delta x_{r,0}$','interpreter','latex','fontsize',15)
xlabel('$|k_1|$','interpreter','latex','fontsize',30)
set(gca,'FontSize',20)
ylim([0 2.9])
