% Function for optimal boundary control
clc; close all; clear variables; 

% DATA
L = 50; % <- Specify domain size
G = 0.5;

%% Space, Time Discretization
N = 501; % number of spatial points
x = linspace(0,L,N); %(x2-x1)/(n-1)
dx = x(2);

T  = 200; % <- Specify final time
dt = 0.5;
tspan = 0:dt:T;

const = [dx N G];

%% Define Desired Traveling Wave Profile hd(t,x) on Domain L, Time T

% Load a pre-defined reference traveling wave (uncomment one of the
% following):

% load('TW1_Ref_L_10_c_0.2.mat'); % defined on domain length 10, wave speed = 0.2
load('TW2_Ref_L_10_c_1.652'); % defined on domain length = 10, wave speed = 1.652
% load('TW3_Ref_L_30_c_2.13'); % defined on domain length = 30, wave speed = 2.13

c = 1.652; % <---- Specify speed of the wave from above 

Nrep = 1; % <--- Specify number of pulses desired on L

hd = zeros(length(tspan),N);
hshift = zeros(1,N);

% Define profile on L (uncomment only if L > 10):
temp = zeros(1,N);
for j = 1:N
    y = rem(x(j),L/Nrep);
    temp(j) = interp1(xref,href,y);
end

for k = 1:length(tperiod)
    for j = 1:N
        y = x(j)-c*tperiod(k);
        if y >= 0
            hshift(j) = interp1(x,temp,y);
        else
            r = rem(abs(y),L);
            hshift(j) = interp1(x,temp,L-r);
        end
    end
    hd(k,:) = hshift;
    clear hshift y r;
end

% Extend hd up to time period T
if T>Tp
    for k = length(tperiod)+1:length(tspan)
        r = mod(k,length(tperiod));
        if r == 0
           hd(k,:) = hd(1,:);
        else
            hd(k,:) = hd(r,:);
        end
    end
end

%% Constant film profile
% If constant film profile is desired instead of traveling wave profile, 
% comment the previous section and uncomment the following lines.

% hd = 1*ones(1,N);
% hd = repmat(hd,length(tspan),1);

%% Initialization

% Initialize Fwd Eq
h0 = 0.5*ones(N-2,1);

% Initiaize Bkwd Eq
p0 = zeros(1,N-4);

% Initialize u(t)
u = 0.5*ones(length(tspan),1);

%% Forward Equation under Initial Control
[~,hfwd] = ode15s(@(t,h)Forward_BC_Func(t,h,const,tspan,u),tspan,h0);

hfwd(:,N-1) = 3/17*(3/2*hfwd(:,N-4)-7*hfwd(:,N-3)+(12-5/6)*hfwd(:,N-2));
hfwd(:,N) = 2/3*(2*hfwd(:,N-1)-1/2*hfwd(:,N-2));

hfwd0 = hfwd;

%% Initial Cost
lambda = 1; % <--- Specify Optimization Parameter

Lh = 0.5*sum(sum((hfwd-hd).^2,2)*dx)*dt;
L0 =  Lh + lambda/2*norm(u)^2*dt;
disp(L0);

%% Backward Equation with initial forward data
[~,p1] = ode15s(@(t,p)Adjoint_BC_Func(T-t,p,const,tspan,hfwd,hd),tspan,p0);

p1_1 = 10/7*p1(:,1)-3/7*p1(:,2); % p1(x=x1,t)
p1_0 = 2/3*(2*p1_1-1/2*p1(:,1)); % p1(x=x0,t)

p2 = flip(p1_0); % p2(t) = p1(0,t)

%% Optimal control computation

% Optimization Parameters
iter = 15; % <------- Specify number of iterations
stepsz = 0.1; % <------ Specify Step size 

Lc = [L0 zeros(1,iter)];

for  i = 1:iter

    % update control
    uc = u - stepsz*(lambda*u + p2);
    
    ind = find(uc<0); % Projected gradient
    if isempty(ind) == 0
        uc(ind)=0;
    end
    
    % Fwd Eq (with updated control)
    [~,hc] = ode15s(@(t,h)Forward_BC_Func(t,h,const,tspan,uc),tspan,h0);

    hc(:,N-1) = 3/17*(3/2*hc(:,N-4)-7*hc(:,N-3)+(12-5/6)*hc(:,N-2));
    hc(:,N) = 2/3*(2*hc(:,N-1)-1/2*hc(:,N-2));

    % New Cost 
    Lh(i) = 0.5*sum(sum((hc(tl:tu,Ll:Lu)-hd(tl:tu,Ll:Lu)).^2,2)*dx)*dt;
    Lc(i) = Lh(i) + lambda/2*norm(uc).^2*dt;
    
    if Lc(i) <  L0
        
        u = uc;
        L0 = Lc(i);
        disp(i); disp(L0);
        hfwd = hc;

        % Bkwd Eq
        [~,p1] = ode15s(@(t,p)Adjoint_BC_Func(T-t,p,const,tspan,hfwd,hd),tspan,p0);

        p1_1 = 10/7*p1(:,1)-3/7*p1(:,2); % p1(x=x1,t)
        p1_0 = 2/3*(2*p1_1-1/2*p1(:,1)); % p1(x=x0,t)
       
        p2 = flip(p1_0); % p2(t) = p1(0,t)

    else
        stepsz = stepsz/2;
    end
end
