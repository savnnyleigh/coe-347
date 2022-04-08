clear all
close all
clc

c=1;        % advective speed
L=4*pi;     % computational domain [0,L]
T=2*2*pi;   % end time
M=0;        % intermediate solutions

fexact='exact.dat';

sigma=0.5; % Courant number
%number of interior points defined in each function

method='forward-upwind';
%method='implicit-central';
%method='beam-warming';
%method='lax-wendroff';

% initial conditions
u0 = @(x) sin(x);  % anonymous function

xx=linspace(0,L,1000);

% solve

out1=wave_solve(c,L,25,sigma,T,M,u0,method);
out2=wave_solve(c,L,50,sigma,T,M,u0,method);
out3=wave_solve(c,L,100,sigma,T,M,u0,method);
out4=wave_solve(c,L,200,sigma,T,M,u0,method);
out5=wave_solve(c,L,400,sigma,T,M,u0,method);

% plot

hold on
grid on
xx=linspace(0,L,1000);

exact(:,2)=u0(xx-out1.TT(2))';
plot(xx,u0(xx-out1.TT(2)),'r-');
plot(out1.x,out1.U(:,2),'--' );
plot(out2.x,out2.U(:,2)) ;
plot(out3.x,out3.U(:,2)),'k-';
plot(out4.x,out4.U(:,2));
plot(out5.x,out5.U(:,2));

axis([0,L,-1.1,1.1]);
xlabel('x');
ylabel('u(x) and numerical solution');
legend('exact','N = 25', 'N = 50', 'N = 100', 'N = 200', 'N=400');
title(sprintf('Time is %f Method is %s Sigma is %f', out1.TT(2), method, sigma));

%{
% dump
fout=sprintf('%s_n%g_sigma%f.dat',method,25,sigma);
dlmwrite(fout,[out1.x',out1.U],'delimiter',' ','precision','%e');
dlmwrite(fexact,[xx',exact],'delimiter',' ','precision','%e');

fout=sprintf('%s_n%g_sigma%f.dat',method,50,sigma);
dlmwrite(fout,[out2.x',out2.U],'delimiter',' ','precision','%e');
dlmwrite(fexact,[xx',exact],'delimiter',' ','precision','%e');

fout=sprintf('%s_n%g_sigma%f.dat',method,100,sigma);
dlmwrite(fout,[out3.x',out3.U],'delimiter',' ','precision','%e');
dlmwrite(fexact,[xx',exact],'delimiter',' ','precision','%e');

fout=sprintf('%s_n%g_sigma%f.dat',method,200,sigma);
dlmwrite(fout,[out4.x',out4.U],'delimiter',' ','precision','%e');
dlmwrite(fexact,[xx',exact],'delimiter',' ','precision','%e');

fout=sprintf('%s_n%g_sigma%f.dat',method,400,sigma);
dlmwrite(fout,[out5.x',out5.U],'delimiter',' ','precision','%e');
dlmwrite(fexact,[xx',exact],'delimiter',' ','precision','%e');

%}