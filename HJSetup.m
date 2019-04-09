function  [x,u0,h,a0,dx,dt,n,Nt] = HJSetup(n,Tf,flag)
% make dx a round number. So want n = 101,
n = n+1;
x = linspace(-3,3,n); dx = x(2)-x(1);  x = x(:); % make it a column vector

h = 0;
a0 = 0.1;
%a0 =0;

% subtract off the small amount of "free" diffusion
a0 = max(0, a0 - dx/2);

% Time step
dt1 = dx;
dt2 = .5*dx^2/a0;
dtmax = min(dt1,dt2);
% round off the time step, making it smaller, as needed
Nt = Tf/dtmax;  Nt = ceil(Nt); dt = Tf/Nt;


if flag == 1
    % function with steps
    u0 = 0*x;
    I1 = abs(x) <= 1; u0(I1) = 2;
    I2 = abs(x) <= .2; u0(I2) = -1;
    %I2 = x <= -.25; u0(I2) = 1;
elseif flag == 2
    u0 = 0*x;
    I1 = abs(x) <= 1; u0(I1) = 2;    
else
    
    % sin function flattened out
    u0 = sin(pi*x); I1 = x >1.5;  u0(I1) = -1; I2 = x < -1.5; u0(I2) = +1;
end