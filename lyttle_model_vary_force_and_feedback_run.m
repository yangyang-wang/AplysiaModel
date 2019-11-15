% Script to run Yangyang's version of the Lyttle model for different forces
% and plot net seaweed intake rate.
%
% To verify results of earlier paper with this model.
%
% PJT w/ HDL 2019-10-24

tic
poolobj = gcp('nocreate');
if max(size(poolobj))==0
    parpool(4) % start a pool of four workers if not already running
end

figure
%num_force=11;
%num_force=401;
%fvec=linspace(0,.05,num_force);
%fvec=logspace(-1,1,num_force);
%fvec=linspace(0,1,num_force);
%fvec=0:.001:.12;
%fvec=0.116:.00005:.117;
fvec=sort([0:.001:.12,0.116:.00005:.117]);
%fvec=0.11620; % near bifurcation point
num_force=length(fvec);
S_rate_in=nan(size(fvec));
parfor j=1:length(fvec)
    disp([j,length(fvec)])
    force=fvec(j);
    M=lyttle_model_vary_force_and_feedback;
    M.tmax=110;
    M.yinit(8)=force; % Force is the 8th "variable".
    M.eps1=1e-3;
    M.eps2=1e-3;
    M.eps3=1e-3;
    M.solve;
    t=M.t;
    S=M.yext(:,7);
    % discard trace for t<=10
    idx_start=find(t>10,1,'first');
    t=t(idx_start:end);
    S=S(idx_start:end);
    
    % find where dS is zero
    dS=diff(S);
    % find last point where dS/dt=0 (this is the open->close transition)
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
        S_rate_in(j)=delta_S/delta_t;
    else
        i1=0; i2=0; delta_t=0; delta_S=0;
        S_rate_in(j)=0;
    end
    %% Plot to check markers are consistent!
    %plot(t,S,'.-')
    %hold on
    %if ~isempty(idx_end_of_flat)
    %    plot(t(i1),S(i1),'ro')
    %    plot(t(i2),S(i2),'go')
    %end
    %beep
    
end
save test_lyttle_vary_force20191115_phase_a.mat
    

%% last figure

beep,beep

figure
plot(fvec,S_rate_in,'+-','LineWidth',3)
set(gca,'FontSize',20)
xlabel('Force Opposing Ingestion')
ylabel('Seaweed/Time')
title('Net Inward Seaweed Rate vs Load')

%% total time
runtime=toc

%save test_lyttle_vary_force20191026b.mat


%delete(poolobj) % optional -- shut down the pool.

% Note for f=0.11620 we have S=0.077817
% but for f=0.11625 we have 0.
