% Use the piecewise iSRC with Piecewise Uniform Rescaling [nu1_close,
% nu1_open] to compute robustness delta_x1/delta_x0 - T1/T0; 
% Then compute its sensitivity to parameter changes

%% vary eps1,eps2,eps3, and fsw, plot the surface of seaweed intake rate -S0. 
% robustness is also calculated, but not ploted

eps1_range = (-6:0.1:-2);
% eps1_range = (-4.8:0.1:-3); % zoom in
fsw_range = (0:0.01:0.1);
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
        
        % Compute Robustness, with perturbation
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

% task_fitness_1storder = dx1_over_dx0 - T1_over_T0; % first order approx of (Seps-S0)/(eps*S0)
Robustness = (Seps-S0)./(eps*S0);  % (Seps-S0)/(eps*S0)

toc 

% save intake-rate.mat eps1_range fsw_range task_fitness_1storder task_fitness T0 Teps S0 Seps delta_x0 delta_xeps

% Robustness/sensitivity to pertbation on fsw
figure
surf(eps1_range,fsw_range,-S0');
shading interp
% hold on; plot3(eps1_range(7),fsw_range(2),task_fitness(7,2),'kp','MarkerEdgeColor','k','MarkerFaceColor',[1,0.5,0.5],'MarkerSize',15)
xlabel('$\rm log10(\varepsilon_0)$','interpreter','latex','fontsize',30)
ylabel('$\rm fsw$','interpreter','latex','fontsize',30,'rot',0)
zlabel('robustness','Interpreter','latex','FontWeight','normal','Fontsize',30)
title('Intake rate $(S_0)$' ,'FontWeight','normal','interpreter','latex','fontsize',25)
set(gca,'FontSize',20)
xlim([eps1_range(1) eps1_range(end)])
ylim([fsw_range(1) fsw_range(end)])

%% robustness to a1
eps2_range = (-6:0.1:-2);
fsw_range = (0:0.01:0.1);
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

% task_fitness_1storder = dx1_over_dx0 - T1_over_T0; % first order approx of (Seps-S0)/(eps*S0)
Robustness = (Seps-S0)./(eps*S0);  % (Seps-S0)/(eps*S0)

toc 

% save intake-rate.mat eps1_range fsw_range task_fitness_1storder task_fitness T0 Teps S0 Seps delta_x0 delta_xeps

% Robustness/sensitivity to pertbation on fsw
figure
surf(eps2_range,fsw_range,-S0');
shading interp
% hold on; plot3(eps1_range(7),fsw_range(2),task_fitness(7,2),'kp','MarkerEdgeColor','k','MarkerFaceColor',[1,0.5,0.5],'MarkerSize',15)
xlabel('$\rm log10(\varepsilon_1)$','interpreter','latex','fontsize',30)
ylabel('$\rm fsw$','interpreter','latex','fontsize',30,'rot',0)
zlabel('robustness','Interpreter','latex','FontWeight','normal','Fontsize',30)
title('Intake rate $(S_0)$' ,'FontWeight','normal','interpreter','latex','fontsize',25)
set(gca,'FontSize',20)
xlim([eps2_range(1) eps2_range(end)])
ylim([fsw_range(1) fsw_range(end)])

%% robustness to a2
% eps3_range = linspace(5e-5,1e-2,20);
eps3_range = (-6:0.1:-2);
fsw_range = (0:0.01:0.1);
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

% task_fitness_1storder = dx1_over_dx0 - T1_over_T0; % first order approx of (Seps-S0)/(eps*S0)
Robustness = (Seps-S0)./(eps*S0);  % (Seps-S0)/(eps*S0)

toc 

% save intake-rate.mat eps1_range fsw_range task_fitness_1storder task_fitness T0 Teps S0 Seps delta_x0 delta_xeps

% Robustness/sensitivity to pertbation on fsw
figure
surf(eps3_range,fsw_range,-S0');
shading interp
% hold on; plot3(eps1_range(7),fsw_range(2),task_fitness(7,2),'kp','MarkerEdgeColor','k','MarkerFaceColor',[1,0.5,0.5],'MarkerSize',15)
xlabel('$\rm log10(\varepsilon_2)$','interpreter','latex','fontsize',30)
ylabel('$\rm fsw$','interpreter','latex','fontsize',30,'rot',0)
zlabel('robustness','Interpreter','latex','FontWeight','normal','Fontsize',30)
title('Intake rate $(S_0)$' ,'FontWeight','normal','interpreter','latex','fontsize',25)
set(gca,'FontSize',20)
xlim([eps3_range(1) eps3_range(end)])
ylim([fsw_range(1) fsw_range(end)])




