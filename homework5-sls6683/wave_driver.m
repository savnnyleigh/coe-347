clear
close all
clc

c=1;        % advective speed
L=2*pi;     % computational domain [0,L]
T=2*2*pi;   % end time
M=0;        % intermediate solutions

fexact='exact.dat';

sigma=1; % Courant number
n=25;       % number of interior points

method1='forward-upwind';
method2='implicit-central';
method3='beam-warming';
method4='lax-wendroff';

% initial conditions
u0 = @(x) sin(x);  % anonymous function

% solve
out1=wave_solve(c,L,n,sigma,T,M,u0,method1);
out2=wave_solve(c,L,n,sigma,T,M,u0,method2);
out3=wave_solve(c,L,n,sigma,T,M,u0,method3);
out4=wave_solve(c,L,n,sigma,T,M,u0,method4);

% plot
hold on
xx=linspace(0,L,1000);
exact(:,2)=u0(xx-out1.TT(2))';
plot(xx,u0(xx-out1.TT(2)),'r-');
plot(out1.x,out1.U(:,2));
plot(out2.x,out2.U(:,1));
plot(out3.x,out3.U(:,2));
plot(out4.x,out4.U(:,2));

axis([0,L,-1.1,1.1]);
xlabel('x');
ylabel('u(x) and numerical solution');
legend('exact','forward-upwind', 'implicit-central', 'beam-warming', 'lax-wendroff');
title(sprintf('Time is %f Sigma is %f',out1.TT(2), sigma));
 
%{
% dump
fout=sprintf('%s_n%g_sigma%f.dat',method,n,sigma);
dlmwrite(fout,[out1.x',out1.U],'delimiter',' ','precision','%e');
dlmwrite(fexact,[xx',exact],'delimiter',' ','precision','%e');

%}