% xdot = [ x_2 ] + [ 0 ] * u
%        [ 0   ]   [ 1 ]
% 
% u(t) \in [-1, 1]
% XT = {(0,0)}
% h = 1, H = 0
% 
% Trajectory starts at (0.3,1), and arrives at (0,0) in minimum time.

clear;
T = 15;
d = 8;
r2 = 5;

% t = sdpvar(1);
% x = sdpvar(2,1,'full');
% u = sdpvar(1);

t = msspoly('t',1);
x = msspoly('x',2);
u = msspoly('u',1);


f = T*[ x(2); 0 ];
g = T*[ 0; 1 ];

% Domains and transitions
y = x;
hX = r2 - y'*y;      
hU = 1 - u^2;       % [-1, 1]
hXT = -y.^2;        % {(0,0)}

% h = 1*x' * x + 20 * u^2;
% H = 0;

h = 1;
H = 0;

x0 = [ 0.3; 1 ];

% Options
options.freeFinalTime = 1;
options.withInputs = 1;
options.svd_eps = 1e4;

[out] = OCPDualSolver(t,x,u,f,g,hX,hU,x0,hXT,h,H,d,options);
% [out] = OCPDualSolverYalmip(t,x,u,f,g,hX,hU,x0,hXT,h,H,d,options);
out.pval*T
