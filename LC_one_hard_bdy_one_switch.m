function [t,y,TE]=LC_one_hard_bdy_one_switch(yinit,tmax,a,b)

%function [t,y,TE]=LC_one_hard_bdy_one_switch(yinit,tmax,a,b)
%
% Matlab function to calculate trajectories of a linear spiral source with
% a hard boundary which forces the system into a stable limit cycle.
% Perhaps the simplest example of how one can get a stable limit cycle with
% a strictly linear vector field, and a sliding boundary.
%
% Use LC_one_hard_bdy_time.m to generate multiple runs and construct a
% "time to liftoff" function (like the asymptotic phase function). 
%
% Here is a way to run it to see a single instance of the limit cycle.
% Change the initial conditions to get different values.  The hard boundary
% is at y=-1, so start at y>=-1 and arbitrary x.
%
%tmax=200;alpha=.2;yinit=[alpha,-1];
%[t,y]=LC_one_hard_bdy(yinit,tmax,alpha);
%figure;plot(y(:,1),y(:,2));grid on;shg
%
% Peter Thomas, CWRU, May 2, 2019 
%
% Added switching surface to create a figure for grant proposal Nov. 21,
% 2019.
%
% Upon crossing x=0 from positive to negative, subtract an amount beta=b>=0 
% from the y-velocity to enact a switching surface (at x=0).

global alpha beta
if nargin < 4, beta=1; else beta=b; end
if nargin < 3, alpha=.2; else alpha=a; end
if nargin < 2, tmax=100; end
if nargin < 1, yinit=[0,-1];end
t=0;
TE=-1; % record exit time (time at which hit liftoff point)
y=yinit;
if (yinit(1)<alpha) && (yinit(2)==-1)
    dom=1; % on wall initially
elseif yinit(1)>0
    dom=0; % right interior domain initially
else
    dom=2; % left interior domain initially
end

%options0=odeset('Events',@dom2_to_wall,'AbsTol',1e-8,'MaxStep',0.01);
%options1=odeset('Events',@wall1_exit,'AbsTol',1e-8,'MaxStep',0.01);
%options2=odeset('Events',@dom0_to_dom2,'AbsTol',1e-8,'MaxStep',0.01);

options0=odeset('Events',@dom0_to_dom2,'AbsTol',1e-8,'MaxStep',0.01);
options1=odeset('Events',@wall1_exit,'AbsTol',1e-8,'MaxStep',0.01);
options2=odeset('Events',@dom2_to_wall,'AbsTol',1e-8,'MaxStep',0.01);

while max(t)<tmax
    switch dom
        case 0 % interior, x>0
            y0=y(end,:)';
            [tnew,ynew,TE,YE,~] = ode45(@dydt0,[max(t),tmax],y0,options0);
            t=[t;tnew];
            y=[y;ynew];
            dom=2; % domain (event index) labels which wall was hit
        case 1 % y=-1 wall
            y0=[y(end,1),-1];
            %y0=[YE(1),-1];
            [tnew,ynew,TE,YE,~] = ode45(@dydt1,[max(t),tmax],y0,options1);
            t=[t;tnew];
            y=[y;ynew];
            dom=0;
            if TE>0
                t=tmax+1; % force end of run upon reaching liftoff point
            end 
        case 2 % interior, x<=0
            y0=y(end,:)';
            [tnew,ynew,TE,YE,~] = ode45(@dydt2,[max(t),tmax],y0,options2);
            t=[t;tnew];
            y=[y;ynew];
            dom=1; 
    end
end

end

function dydt=dydt0(~,y) % flow on the interior, x>0
    global alpha
    dydt=[alpha,-1;1,alpha]*y;
end

function dydt=dydt1(~,y) % flow on the wall y=-1
    global alpha
    dydt=[alpha,-1;1,alpha]*y;
    dydt(2)=max(dydt(2),0); % only allow positive dy/dt or else zero
end

function dydt=dydt2(~,y) % flow on the interior, x<=0
    global alpha beta
    dydt=[alpha,-1;1,alpha]*y - beta*[0;1];
end

function [value,isterminal,direction]=dom2_to_wall(~,y)
    value=y(2)+1;    % when y crosses -1 from above (wall)
    isterminal=1; % stop integration and return
    direction=-1; % "value" should be decreasing
end

function [value,isterminal,direction]=wall1_exit(~,y)
    % when the *unconstrained* value of dy/dt increases through zero, return
    global alpha
    dydt=[alpha,-1;1,alpha]*y; 
    value=dydt(2);
    isterminal=1;
    direction=1;
end

function [value,isterminal,direction]=dom0_to_dom2(~,y)
    value=-y(1); % when x crosses 0 from above
    isterminal=1; % stop integration and return
    direction=1; % "value" should be decreasing
end
