% compute intake rate, period, net grasper change, and robustness to
% perturbation on fsw: fsw -> fsw+eps as we vary fsw (load size) and eps1
% (sensory feedback strength to a2 neural pool)

% delta_x1/delta_x0-T1/T0
%        where dx0 is the net change in grasper position xr during closing phase under fsw
%              delta_xeps is .... under fsw+eps, dxeps=dx0+eps*dx1

% % in this range: model is at HC mode
% eps3_range = (4e-5:1e-5:2e-4);
% fsw_range = (0.005:0.01:0.1);

% bigger range but sparser points 
% eps3_range = [(1e-5:1e-5:1e-4) linspace(1e-4,1e-3,10)]; 
% eps3_range = [(1e-7:1e-6:2e-5) linspace(2e-5,1e-3,10)];
% fsw_range = (0:0.02:0.5);

% % eps3_range = (1e-6:2e-6:6e-5);
% eps3_range = (6e-5:4e-5:8e-4);
% fsw_range = 0.045;

eps3_range = linspace(5e-5,1e-3,10);
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
        
        eps3=eps3_range(i);
        fsw = fsw_range(j);
        
        % no perturbation, with parameters eps3 and fsw
        model1 = lyttle_model('tmax',20,'eps3',eps3); % assign new eps3 value
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

        if isempty(model_pert.t_open_to_close)
            continue
        end
        
        Teps(i,j)= model_pert.t_open_to_close(1);
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

task_fitness_1storder = dx1_over_dx0 - T1_over_T0; % first order approx of (Seps-S0)/(eps*S0)
task_fitness = (Seps-S0)./(eps*S0);  % (Seps-S0)/(eps*S0)

toc 

% save intake-rate.mat eps3_range fsw_range task_fitness_1storder task_fitness T0 Teps S0 Seps delta_x0 delta_xeps

%% Robustness/sensitivity to pertbation on fsw
figure
surf(eps3_range,fsw_range,task_fitness');
shading interp
% hold on; plot3(eps3_range(7),fsw_range(2),task_fitness(7,2),'kp','MarkerEdgeColor','k','MarkerFaceColor',[1,0.5,0.5],'MarkerSize',15)
xlabel('$\varepsilon_2$','interpreter','latex','fontsize',30)
ylabel('$\rm fsw$','interpreter','latex','fontsize',30,'rot',0)
zlabel('robustness','Interpreter','latex','FontWeight','normal','Fontsize',30)
title('robustness to increase in fsw' ,'FontWeight','normal','interpreter','latex','fontsize',25)
set(gca,'FontSize',20)
xlim([eps3_range(1) eps3_range(end)])
ylim([fsw_range(1) fsw_range(end)])

view([0 90]);
colorbar

% seewead intake rate w.r.t. eps3 and fsw
figure
surf(eps3_range,fsw_range,-S0'); % need to add the negative sign! because:  -delta_x0 = seaweed intake
shading interp
% hold on; plot3(eps3_range(7),fsw_range(2),-S0(7,2),'kp','MarkerEdgeColor','k','MarkerFaceColor',[1,0.5,0.5],'MarkerSize',15)
xlabel('$\varepsilon_2$','interpreter','latex','fontsize',30)
ylabel('$\rm fsw$','interpreter','latex','fontsize',30,'rot',0)
zlabel('intake rate','Interpreter','latex','FontWeight','normal','Fontsize',30)
title('seaweed intake rate' ,'FontWeight','normal','interpreter','latex','fontsize',25)
set(gca,'FontSize',20)
xlim([eps3_range(1) eps3_range(end)])
ylim([fsw_range(1) fsw_range(end)])

view([0 90]);
colorbar

%% 3d
figure
subplot(2,2,1)
surf(eps3_range,fsw_range,-S0'); % need to add the negative sign! because:  -delta_x0 = seaweed intake
shading interp
% hold on; plot3(eps3_range(7),fsw_range(2),-S0(7,2),'kp','MarkerEdgeColor','k','MarkerFaceColor',[1,0.5,0.5],'MarkerSize',15)
xlabel('$\varepsilon_2$','interpreter','latex','fontsize',30)
ylabel('$\rm fsw$','interpreter','latex','fontsize',30,'rot',0)
title('seaweed intake rate $-\Delta x_r/T_0$' ,'FontWeight','normal','interpreter','latex','fontsize',25)
set(gca,'FontSize',20)
xlim([eps3_range(1) eps3_range(end)])
ylim([fsw_range(1) fsw_range(end)])


subplot(2,2,2)
surf(eps3_range,fsw_range,task_fitness');
shading interp
% hold on; plot3(eps3_range(7),fsw_range(2),task_fitness(7,2),'kp','MarkerEdgeColor','k','MarkerFaceColor',[1,0.5,0.5],'MarkerSize',15)
xlabel('$\varepsilon_2$','interpreter','latex','fontsize',30)
ylabel('$\rm fsw$','interpreter','latex','fontsize',30,'rot',0)
title('robustness to increase in fsw' ,'FontWeight','normal','interpreter','latex','fontsize',25)
set(gca,'FontSize',20)
xlim([eps3_range(1) eps3_range(end)])
ylim([fsw_range(1) fsw_range(end)])


subplot(2,2,3)
surf(eps3_range,fsw_range,T0'); % need to add the negative sign! because:  -delta_x0 = seaweed intake
shading interp
% hold on; plot3(eps3_range(7),fsw_range(2),T0(7,2),'kp','MarkerEdgeColor','k','MarkerFaceColor',[1,0.5,0.5],'MarkerSize',15)
xlabel('$\varepsilon_2$','interpreter','latex','fontsize',30)
ylabel('$\rm fsw$','interpreter','latex','fontsize',30,'rot',0)
title('period $T_0$' ,'FontWeight','normal','interpreter','latex','fontsize',25)
set(gca,'FontSize',20)
xlim([eps3_range(1) eps3_range(end)])
ylim([fsw_range(1) fsw_range(end)])


subplot(2,2,4)
surf(eps3_range,fsw_range,delta_x0'); % need to add the negative sign! because:  -delta_x0 = seaweed intake
shading interp
% hold on; plot3(eps3_range(7),fsw_range(2),delta_x0(7,2),'kp','MarkerEdgeColor','k','MarkerFaceColor',[1,0.5,0.5],'MarkerSize',15)
xlabel('$\varepsilon_2$','interpreter','latex','fontsize',30)
ylabel('$\rm fsw$','interpreter','latex','fontsize',30,'rot',0)
title('net change in $x_r: \Delta x_r$' ,'FontWeight','normal','interpreter','latex','fontsize',25)
set(gca,'FontSize',20)
xlim([eps3_range(1) eps3_range(end)])
ylim([fsw_range(1) fsw_range(end)])


%% 2d
figure
subplot(2,2,1)
surf(eps3_range,fsw_range,-S0'); % need to add the negative sign! because:  -delta_x0 = seaweed intake
shading interp
% hold on; plot3(eps3_range(7),fsw_range(2),-S0(7,2),'kp','MarkerEdgeColor','k','MarkerFaceColor',[1,0.5,0.5],'MarkerSize',15)
xlabel('$\varepsilon_2$','interpreter','latex','fontsize',30)
ylabel('$\rm fsw$','interpreter','latex','fontsize',30,'rot',0)
title('seaweed intake rate $-\Delta x_r/T_0$' ,'FontWeight','normal','interpreter','latex','fontsize',25)
set(gca,'FontSize',20)
xlim([eps3_range(1) eps3_range(end)])
ylim([fsw_range(1) fsw_range(end)])
view([0 90]);
colorbar

subplot(2,2,2)
surf(eps3_range,fsw_range,task_fitness');
shading interp
% hold on; plot3(eps3_range(7),fsw_range(2),task_fitness(7,2),'kp','MarkerEdgeColor','k','MarkerFaceColor',[1,0.5,0.5],'MarkerSize',15)
xlabel('$\varepsilon_2$','interpreter','latex','fontsize',30)
ylabel('$\rm fsw$','interpreter','latex','fontsize',30,'rot',0)
title('robustness to increase in fsw' ,'FontWeight','normal','interpreter','latex','fontsize',25)
set(gca,'FontSize',20)
xlim([eps3_range(1) eps3_range(end)])
ylim([fsw_range(1) fsw_range(end)])
view([0 90]);
colorbar

subplot(2,2,3)
surf(eps3_range,fsw_range,T0'); % need to add the negative sign! because:  -delta_x0 = seaweed intake
shading interp
% hold on; plot3(eps3_range(7),fsw_range(2),T0(7,2),'kp','MarkerEdgeColor','k','MarkerFaceColor',[1,0.5,0.5],'MarkerSize',15)
xlabel('$\varepsilon_2$','interpreter','latex','fontsize',30)
ylabel('$\rm fsw$','interpreter','latex','fontsize',30,'rot',0)
title('period $T_0$' ,'FontWeight','normal','interpreter','latex','fontsize',25)
set(gca,'FontSize',20)
xlim([eps3_range(1) eps3_range(end)])
ylim([fsw_range(1) fsw_range(end)])
view([0 90]);
colorbar

subplot(2,2,4)
surf(eps3_range,fsw_range,delta_x0'); % need to add the negative sign! because:  -delta_x0 = seaweed intake
shading interp
% hold on; plot3(eps3_range(7),fsw_range(2),delta_x0(7,2),'kp','MarkerEdgeColor','k','MarkerFaceColor',[1,0.5,0.5],'MarkerSize',15)
xlabel('$\varepsilon_2$','interpreter','latex','fontsize',30)
ylabel('$\rm fsw$','interpreter','latex','fontsize',30,'rot',0)
title('net change in $x_r: \Delta x_r$' ,'FontWeight','normal','interpreter','latex','fontsize',25)
set(gca,'FontSize',20)
xlim([eps3_range(1) eps3_range(end)])
ylim([fsw_range(1) fsw_range(end)])
view([0 90]);
colorbar