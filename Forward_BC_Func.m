% Function for thin-film equation with boundary control 
% Central-difference scheme with no ghost points
% Updated: 3/8/22

function dhdt = Forward_BC_Func(t,h0,C,tspan,u0)

dx = C(1);
N = C(2);
G = C(3);

q_0 = interp1(tspan,u0,t);

dhdt = zeros(N-2,1);
q = zeros(N-1,1);

h = zeros(N,1);
h(1:N-2) = h0; % DE solved for grid points 3 to N-2

% Boundary points (using midpoint values at previous time step)
h(N-1) = 3/17*(3/2*h(N-4)-7*h(N-3)+(12-5/6)*h(N-2));
h(N) = 2/3*(2*h(N-1)-1/2*h(N-2));

% i = 2
q(1) = q_0;

q(3) = h(3)^3*(G + (h(4)-h(2))/(2*dx) + (h(5)-2*h(4)+2*h(2)-h(1))/(2*dx^3));

dhdt(2) =  -(q(3)-q(1))/(2*dx);

% i = 3
hxxx = (-3*h(6)+14*h(5)-24*h(4)+18*h(3)-5*h(2))/(2*dx^3); % fwd diff
q(2) = h(2)^3*(G + (h(3)-h(1))/(2*dx) + hxxx);

q(4) = h(4)^3*(G + (h(5)-h(3))/(2*dx) + (h(6)-2*h(5)+2*h(3)-h(2))/(2*dx^3));

dhdt(3) =  -(q(4)-q(2))/(2*dx);

for i = 4:N-3

    q(i+1) = h(i+1)^3*(G + (h(i+2)-h(i))/(2*dx) + ...
        (h(i+3)-2*h(i+2)+2*h(i)-h(i-1))/(2*dx^3));

    dhdt(i) =  -(q(i+1)-q(i-1))/(2*dx);

end

% i = N-2
hxxx = (5/2*h(N-1)-9*h(N-2)+12*h(N-3)-7*h(N-4)+3/2*h(N-5))/(dx^3); 

q(N-1) = h(N-1)^3*(G + (h(N)-h(N-2))/(2*dx) + hxxx);

dhdt(N-2) = -(q(N-1)-q(N-3))/(2*dx);

end