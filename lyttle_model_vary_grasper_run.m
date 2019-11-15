% Script to run Yangyang's version of the Lyttle model for different
% thresohld geometries and plot net seaweed intake rate.
%
% PJT w/ help from YYW 2019-11-01

[status,hostname]=system('echo $HOSTNAME'); % are we on phase or laptop?

tic
poolobj = gcp('nocreate');
if max(size(poolobj))==0
    if strcmp(hostname(1:end-1),'phase.MATH.CWRU.Edu')
        parpool(6) % start a pool of six workers if not already running
    else
        parpool(4) % start a pool of four workers if not already running
    end
end

figure
force=0.1; % default = force=.01;
nthresh=75;
clist=linspace(-.5,1.5,nthresh); % list of thresholds
nangle=60;
alist=linspace(0,2*pi,nangle); % list of angles
[c,angle]=meshgrid(clist,alist);
S_rate_in=nan(nangle,nthresh);

% If restarting partway through, load partial results here:
%load lyttle_model_vary_grasper20191106a

forcelist=[.02,.03,.04,.05,.06,.07,.08,.09,.11]; % already did 0.01 and 0.10
for force=forcelist
    for i=1:nthresh % 1:nthresh
        parfor j=1:nangle
            disp([j,i,nangle,nthresh])
            M=lyttle_model_vary_grasper;
            M.tmax=40;
            M.yinit(8)=force; % Force is the 8th "variable".
            M.ocangle=alist(j);
            M.octhresh=clist(i)/sqrt(2);
            M.solve;
            t=M.t;
            S=M.yext(:,7);
            % discard trace for t<=10
            idx_start=find(t>10,1,'first');
            t=t(idx_start:end);
            S=S(idx_start:end);
            
            % Find where dS is zero.
            dS=diff(S);
            % Find last point where dS/dt=0 (this is the open->close transition).
            idx_end_of_flat=2+find(...
                (dS(1:end-2)==0).*...
                (dS(2:end-1)==0).*...
                (dS(3:end)>1e-10));
            if ~isempty(idx_end_of_flat)
                i1=idx_end_of_flat(1);
                i2=idx_end_of_flat(end);
                % Net time
                delta_t=t(i2)-t(i1);
                % Net seaweed movement (let positive be inwards)
                delta_S=S(i1)-S(i2);
                % Rate
                S_rate_in(j,i)=delta_S/delta_t;
            else
                i1=0; i2=0; delta_t=0; delta_S=0;
                S_rate_in(j,i)=0;
            end
            
        end
        save(['lyttle_model_vary_grasper20191108_force_',num2str(force)])
    end
end

%% last figure

figure
pcolor(c,angle/pi,S_rate_in)
xlabel('Threshold')
ylabel('Angle/\pi')

%pcolor(angle'/pi,c',S_rate_in') % transpose axes
%ylabel('Threshold')
%xlabel('Angle/\pi')

set(gca,'FontSize',20)
title('Net Inward Seaweed Rate vs Grasper Geometry')
colormap gray
colorbar

%% total time
runtime=toc;
disp(runtime)

%delete(poolobj) % optional -- shut down the pool.

