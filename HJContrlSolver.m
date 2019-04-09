% code for basic finite differences for
% Hamilton-Jacobi solver for
% eikonal with diffusion equation
% The equation is 
% u_t = max(|ux| -1, 0) + a_0 u_xx
% Here running cost is 1.  Can move at speed 1.  Want to maximize, payoff
% at final time is u_0.  Also have a (possibly zero) small amount of
% randomness in the u_xx term.

% Set up the domain
n = 200;
Tf = .5;
flag = 2;
% Set up the initial data
[x,u0,h,a0,dx,dt,n,Nt] = HJSetup(n,Tf,flag);


%% Now iterate to solve
u = u0;

% Trick to keep number of plots small:
nplots = 6;
aa = max(1,floor(Nt/nplots));
mu = dt/dx;

for jj = 1: Nt
    [Hux,H2ux, uxx2] = HJFD(u,dx);
    u = (1-mu)*u + mu*H2ux + dt*a0*uxx2;
%    if mod(jj,aa) == 0 % whether to plot, only want to plot a few times
%    figure(2), plot(x,u,'-*',x,u0); pause(0.5)    
%    end
end
figure(2), plot(x,u,'-*',x,u0);

u = u0;
for jj = 1: Nt
    [Hux,H2ux, uxx2] = HJFD(u,dx);  
    u = (1-mu)*u + mu*Hux + dt*a0*uxx2;
%    if mod(jj,aa) == 0 % whether to plot, only want to plot a few times
%    figure(3), plot(x,u,'-*',x,u0); pause(0.5)    
%    end
end
figure(3), plot(x,u,'-*',x,u0);

u = u0;
for jj = 1: Nt
    [Hux,H2ux, uxx2] = HJFD(u,dx);  
    u = (1-mu)*u + mu*H2ux;
%    if mod(jj,aa) == 0 % whether to plot, only want to plot a few times
%    figure(4), plot(x,u,'-*',x,u0); pause(0.5)    
%    end
end
figure(4), plot(x,u,'-*',x,u0);

u = u0;
for jj = 1: Nt
    [Hux,H2ux, uxx2] = HJFD(u,dx);  
    u = (1-mu)*u + mu*Hux;
%    if mod(jj,aa) == 0 % whether to plot, only want to plot a few times
%    figure(5), plot(x,u,'-*',x,u0); pause(0.5)    
%    end
end
figure(5), plot(x,u,'-*',x,u0);

