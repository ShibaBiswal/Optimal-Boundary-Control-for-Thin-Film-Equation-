% Function for the adjoint equation simulation
% Central-difference scheme with no ghost points
% Updated: 11/8/22

function dpdt = Adjoint_BC_Func(t,q,C,tfwd,hfwd,hdr)

dx = C(1);
N = C(2);
G = C(3);

dpdt = zeros(N-4,1);

h = interp1(tfwd,hfwd,t);
hd = interp1(tfwd,hdr,t);

p = zeros(N,1);
p(3:N-2) = q; % DE solved for grid points 3 to N-2

% Boundary points (using midpoint values at previous time step)
p(2) = 10/7*p(3)-3/7*p(4);
p(1) = 2/3*(2*p(2)-1/2*p(3));

p(N-1) = (3/2*p(N-4)-7*p(N-3)+(12-1/3*(5/2-dx^3))*p(N-2))/(-9+4/3*(5/2-dx^3));
p(N) = 1/3*(4*p(N-1)-p(N-2));

%--------------------------------------------------------------------------
i = 3;
dhdx = (h(i+1)-h(i-1))/(2*dx);
d3hdx3 = (h(i+2)-2*h(i+1)+2*h(i-1)-h(i-2))/(2*dx^3);

dpdx = (p(i+1)-p(i-1))/(2*dx);               % i
dpdx_p = (p(i+2)-p(i))/(2*dx);               % i+1
dpdx_pp = (p(i+3)-p(i+1))/(2*dx);            % i+2 
dpdx_n = (p(i)-p(i-2))/(2*dx);               % i-1
dpdx_nn = (-3*p(i-2)+4*p(i-1)-p(i))/(2*dx); % Fwd diff centered at i-2

dpdt(i-2) = h(i)-hd(i) + dpdx*( 3*h(i)^2 * ( G + dhdx + d3hdx3 ) ) ...
    - ( dpdx_p*h(i+1)^3 - dpdx_n*h(i-1)^3 )/(2*dx) ...
    - ( dpdx_pp*h(i+2)^3 - 2*dpdx_p*h(i+1)^3 ...
        + 2*dpdx_n*h(i-1)^3 - dpdx_nn*h(i-2)^3 )/(2*dx^3);

%--------------------------------------------------------------------------
for  i = 4:N-3

    dhdx = (h(i+1)-h(i-1))/(2*dx);
    d3hdx3 = (h(i+2)-2*h(i+1)+2*h(i-1)-h(i-2))/(2*dx^3);

    dpdx = (p(i+1)-p(i-1))/(2*dx);      % i
    dpdx_p = (p(i+2)-p(i))/(2*dx);      % i+1
    dpdx_pp = (p(i+3)-p(i+1))/(2*dx);   % i+2 
    dpdx_n = (p(i)-p(i-2))/(2*dx);      % i-1
    dpdx_nn = (p(i-1)-p(i-3))/(2*dx);   % i-2

    dpdt(i-2) = h(i)-hd(i) + dpdx *( 3*h(i)^2 * ( G + dhdx + d3hdx3 ) ) ...
        - ( dpdx_p*h(i+1)^3 - dpdx_n*h(i-1)^3 )/(2*dx) ...
        - ( dpdx_pp*h(i+2)^3 - 2*dpdx_p*h(i+1)^3 ...
            + 2*dpdx_n*h(i-1)^3 - dpdx_nn*h(i-2)^3 )/(2*dx^3); 
end

%--------------------------------------------------------------------------
i = N-2;

dhdx = (h(i+1)-h(i-1))/(2*dx);
d3hdx3 = (h(i+2)-2*h(i+1)+2*h(i-1)-h(i-2))/(2*dx^3);
    
dpdx = (p(i+1)-p(i-1))/(2*dx);               % i
dpdx_p = (p(i+2)-p(i))/(2*dx);               % i+1
dpdx_pp = (p(i)-4*p(i+1)+3*p(i+2))/(2*dx);  % Bkwrd diff centered at i+2 
dpdx_n = (p(i)-p(i-2))/(2*dx);               % i-1
dpdx_nn = (p(i-1)-p(i-3))/(2*dx);            % i-2

dpdt(i-2) = h(i)-hd(i) + dpdx *( 3*h(i)^2 * ( G + dhdx + d3hdx3 ) ) ...
    - ( dpdx_p*h(i+1)^3 - dpdx_n*h(i-1)^3 )/(2*dx) ...
    - ( dpdx_pp*h(i+2)^3 - 2*dpdx_p*h(i+1)^3 ...
        + 2*dpdx_n*h(i-1)^3 - dpdx_nn*h(i-2)^3 )/(2*dx^3);

end