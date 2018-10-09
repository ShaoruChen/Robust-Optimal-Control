function [out] = OCPDualSolver(t,x,u,f,g,hX,hU,x0,hXT,h,H,d,options)
%OCPDUALSOLVER Summary of this function goes here
%   Detailed explanation goes here

if mod(d,2) ~= 0
    warning('d is not even. Using d+1 instead.');
    d = d+1;
end


max_m = 0;

m = length(u);
if m > max_m
    max_m = m;
end
if (length(f) ~= length(x)) || (size(g,1) ~= length(x))
    error('Inconsistent matrix size.');
end
if size(g,2) ~= m
    error('Inconsistent matrix size.');
end


svd_eps = 1e3;
if ~exist('options','var'), options = struct(); end
if isfield(options, 'svd_eps'), svd_eps = options.svd_eps; end
if ~isfield(options, 'freeFinalTime'), options.freeFinalTime = 0; end
if ~isfield(options, 'withInputs'), options.withInputs = 0; end


T = 1;
hT = t*(T-t);

prog = spotsosprog;
prog = prog.withIndeterminate( t );
prog = prog.withIndeterminate( x );
prog = prog.withIndeterminate( u );

% create v(i)
vmonom = monomials( [ t; x], 0:d );
[ prog, v, ~ ] = prog.newFreePoly( vmonom );

    
vT = subs( v, t, T );
dvdt = diff( v, t );
dvdx = diff( v, x );
Lfv = dvdt + dvdx*f;
Lgv= dvdx*g;
Lv = Lfv + Lgv*u;

obj = 0;

% Lv_i + h_i >= 0                   Dual: mu
prog = sosOnK( prog, Lv + h, ...
               [ t; x; u ], [ hT; hX; hU ], d);
mu_idx = size( prog.sosExpr, 1 );

% v(T,x) <= H_i(x)                  Dual: muT
if ~isempty( hXT )
    if options.freeFinalTime
        prog = sosOnK( prog, H - v, [ t; x ], [ hT; hXT ], d );
    else
        prog = sosOnK( prog, H - vT, x, hXT, d );
    end
end
    
    
if ~isempty( x0 )
    obj = obj + subs( v, [ t; x ], [ 0; x0 ] );
end

% set options
spot_options = spot_sdp_default_options();
spot_options.verbose = 1;
if isfield(options, 'solver_options')
    spot_options.solver_options = options.solver_options;
end

%solve
tic;
[sol, y, dual_basis] = prog.minimize( -obj, @spot_mosek, spot_options );

out.time = toc;

out.pval = double(sol.eval(obj));
out.sol = sol;

end

