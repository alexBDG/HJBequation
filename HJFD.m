function [Hux,H2ux, uxx2] = HJFD(u,dx)
% the Hamiltonian is 
%  H1(p) = |p| -1
%  H2(p) = max(0,|p|-1) = max(0,H1(p))
% or
%
% shift the vectors using the appropriate BC
[uB, uF] = shiftNeumann(u);
% 
%Hux = max(max(uB-dx,uF-dx),u);
%H2ux = max(u,Hux);
Hux = max(uB, uF);
H2ux = max(max(uB-dx,uF-dx), u);
uxx2 = .5*(uF + uB -2*u)/dx^2;

function [uB, uF] = shiftNeumann(u)
u = u(:); n = length(u);
% duplicating the first term leads to u_x = 0 at left boundary
uB = [u(1); u(1:n-1)];
uF = [u(2:n); u(n)];

function [uB, uF] = shiftDirichlet(u)
u = u(:); n = length(u);
% zero Dirichlet BC at endpoints
uB = [0; u(1:n-1)];
uF = [u(2:n);0];