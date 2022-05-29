function [t_isrc,x_isrc]=lyttle_iSRC(pt,ICisClose,isUniform)
% This function computes and plots ISRC curve

if nargin<1, ICisClose=true;end
if nargin<2, pt=false; end
if nargin<3, isUniform=false;end


model0 = lyttle_model; model0.tmax=20;model0.solve;


eps = 0.02; % fsw -> fsw + eps

% get the IC depending on starting from close or open
if ICisClose
    
    % restart at close: y_open_to_close
    xinit = model0.y_open_to_close(end,1:8);
    
    fsw = model0.xinit(8);
    model = lyttle_model('xinit',xinit,'tmax',10); % start at y_open_to_close
    model.solve
    T0=model.t_open_to_close(1);
    T0_close = model.t_close_to_open(1); % Compute the total time spent during closing
    T0_open=T0-T0_close;
    
    %  Find IC for iSRC
%     eps = 0.02; % fsw -> fsw + eps
    % we are free to choose the difference at the open_to_close point to be the IC for var equation with perturbation on fsw
    model0.xinit(8)=fsw+eps;model0.tmax=10;
    model0.solve;
    xinit_pert = model0.y_open_to_close(end,1:8);   % this is the initial cond for perturbed LC
    % Initial value of iSRC is the difference between open_to_close points
    vinit = (xinit_pert(1:6)-xinit(1:6))/eps;
    
    
elseif ~ICisClose
    % restart at y_close_to_open
    xinit = model0.y_close_to_open(end,1:8);
    
    fsw = model0.xinit(8);
    model = lyttle_model('xinit',xinit,'tmax',10); % start at y_open_to_close
    model.solve
    T0=model.t_close_to_open(1);
    
    T0_open = model.t_open_to_close(1); % Compute the total time spent during open
    T0_close = T0 - T0_open;             % total time spent during close
    
    
    %  Find IC for iSRC
%     eps = 0.02; % fsw -> fsw + eps
    % we are free to choose the difference at the open_to_close point to be the IC for var equation with perturbation on fsw
    model0.xinit(8)=fsw+eps;model0.tmax=10;
    model0.solve;
    xinit_pert = model0.y_close_to_open(end,1:8);   % this is the initial cond for perturbed LC
    % Initial value of iSRC is the difference between open_to_close points
    vinit = (xinit_pert(1:6)-xinit(1:6))/eps;
    
end

% % T1=8.0778; % computed from T1=lyttle_iPRC; this is the linear shift in full period, should approach T1 = (Teps-T0)/eps as eps->0
T1=8.1225; % from numeric

t_isrc=2*T0;

% create the model depending on the way of rescaling and the starting point
if isUniform
    
    % Find the global nu1's for iSRC
    nu1 = T1/T0;
    
    % Compute the unperturbed solution and the iSRC with piecewise nu's found from local timing response curve
    model = lyttle_model( 'xinit', xinit, 'vinit', vinit, ...
        'tmax', t_isrc, 'nu', [nu1 nu1]);
    
elseif ~isUniform
    
    % nu1's for piecewise iSRC
%     T1_close=5.1817; % computed from T1_close = lyttle_ltrc;
    T1_close=5.1848; % from numeric
    nu1_close=T1_close/T0_close;      % nu in interior
    nu1_open=(T1-T1_close)/T0_open; % compute the nu in boundary
    
    % Compute the unperturbed solution and the iSRC with piecewise nu's found from local timing response curve
    model = lyttle_model( 'xinit', xinit, 'vinit', vinit, ...
        'tmax', t_isrc, 'nu', [nu1_close,nu1_open]);
    
end

model.solve

t_isrc = model.t;
x_isrc = model.yext(:,9:14);

if pt
    model.plot_var_horizontal;
end
end