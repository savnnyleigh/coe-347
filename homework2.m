%% Exact Solution
syms y(x)
ode = diff(y,x) == -50*(y-cos(x));
cond = y(0) == 0;
ySol(x) = dsolve(ode,cond);

%% Explicit Euler (forward)


[x1,y1] = feuler(0, 0.033,30);
[x2,y2] = feuler(0, 0.04,25);
[x3,y3] = feuler(0, 0.1,10);
[x4,y4] = feuler(0, 0.02,50);

figure
hold on
grid on

x = 0:.0001:1;
exact = (50*2501^(1/2)*cos(- atan(1/50) + x))/2501 - (2500*exp(-50.*x))/2501;

plot(x3,y3, 'Color',[0 0.4470 0.7410]);
plot(x2,y2, 'Color', [0.8500 0.3250 0.0980]);
plot(x1,y1,'Color', [0.4660 0.6740 0.1880]);
plot(x4,y4,'Color', [0.9290 0.6940 0.1250]);

plot(x,exact, 'Color', [0.4940 0.1840 0.5560]);

legend('h = 0.1, 10 steps','h = 0.04, 25 steps', 'h = 0.033, 30 steps', 'h = 0.02, 50 steps', 'Exact Solution');
xlim([0 1]);
ylim([0 1.8]);
xlabel('x');
ylabel('y');
title('Explicit Euler Solutions');

hold off

e = zeros(1,4);


%% Implicit Euler (backwards)
[x1,y1] = beuler(0, 0.033,30);
[x2,y2] = beuler(0, 0.1,10);
[x3,y3] = beuler(0,.5 ,2);
[x4,y4] = beuler(0, 0.02,50);

figure
hold on
grid on

x = 0:.0001:1;
exact = (50*2501^(1/2)*cos(- atan(1/50) + x))/2501 - (2500*exp(-50.*x))/2501;

plot(x3,y3, 'Color',[0 0.4470 0.7410]);
plot(x2,y2, 'Color', [0.8500 0.3250 0.0980]);
plot(x1,y1,'Color', [0.4660 0.6740 0.1880]);
plot(x4,y4,'Color', [0.9290 0.6940 0.1250]);

plot(x,exact, 'Color', [0.4940 0.1840 0.5560]);

legend('h = 0.5, 2 steps','h = 0.1, 10 steps', 'h = 0.033, 30 steps', 'h = 0.02, 50 steps', 'Exact Solution');
xlim([0 1]);
ylim([0 1.8]);
xlabel('x');
ylabel('y');
title('Implicit Euler Solutions');

hold off

%% Midpoint

[x1,y1] = midpoint(0,0.033,30);
[x2,y2] = midpoint(0, 0.04,25);
[x3,y3] = midpoint(0, 0.1,10);
[x4,y4] = midpoint(0, 0.02,50);

figure
hold on
grid on

x = 0:0.0001:1;
exact = (50*2501^(1/2)*cos(- atan(1/50) + x))/2501 - (2500*exp(-50.*x))/2501;

%plot(x3,y3, 'Color',[0 0.4470 0.7410]);
%plot(x2,y2, 'Color', [0.8500 0.3250 0.0980]);
plot(x1,y1,'Color', [0.4660 0.6740 0.1880]);
plot(x4,y4,'Color', [0.9290 0.6940 0.1250]);

plot(x,exact, 'Color', [0.4940 0.1840 0.5560]);

legend('h = 0.1, 10 steps','h = 0.04, 25 steps', 'h = 0.033, 30 steps', 'h = 0.02, 50 steps', 'Exact Solution');
xlim([0 1]);
ylim([0 1.8]);
xlabel('x');
ylabel('y');
title('Midpoint Solutions');

hold off

 
%% Trapezoidal


[x1,y1] = trapezoid(0,0.033,30);
[x2,y2] = trapezoid(0, 0.04,25);
[x3,y3] = trapezoid(0, 0.1,10);
[x4,y4] = trapezoid(0, 0.02,50);

figure('Name', 'Trapezoidal')
hold on
grid on

x = 0:0.0001:1;
exact = (50*2501^(1/2)*cos(- atan(1/50) + x))/2501 - (2500*exp(-50.*x))/2501;

plot(x3,y3, 'Color',[0 0.4470 0.7410]);
plot(x2,y2, 'Color', [0.8500 0.3250 0.0980]);
plot(x1,y1,'Color', [0.4660 0.6740 0.1880]);
plot(x4,y4,'Color', [0.9290 0.6940 0.1250]);

plot(x,exact, 'Color', [0.4940 0.1840 0.5560]);

legend('h = 0.1, 10 steps','h = 0.04, 25 steps', 'h = 0.033, 30 steps', 'h = 0.02, 50 steps', 'Exact Solution');
xlim([0 1]);
xlabel('x');
ylabel('y');
title('Trapezoidal Solutions');

hold off


%% Adams-Bashforth 2

[x1,y1] = AB2(0, 0.033,30);
[x2,y2] = AB2(0, 0.025,40);
[x3,y3] = AB2(0, 0.01,100);
[x4,y4] = AB2(0, 0.02,50);

figure('Name', 'AB2')
hold on
grid on

x = 0:0.0001:1;
exact = (50*2501^(1/2)*cos(- atan(1/50) + x))/2501 - (2500*exp(-50.*x))/2501;


plot(x1,y1,'Color', [0.4660 0.6740 0.1880]);
plot(x2,y2, 'Color', [0.8500 0.3250 0.0980]);
plot(x4,y4,'Color', [0.9290 0.6940 0.1250]);
plot(x3,y3, 'Color',[0 0.4470 0.7410]);

plot(x,exact, 'Color', [0.4940 0.1840 0.5560]);

legend('h = 0.033, 30 steps', 'h = 0.025, 40 steps', 'h = 0.02, 50 steps', 'h = 0.01, 100 steps', 'Exact Solution');
xlim([0 1]);
ylim([0 1.8]);
xlabel('x');
ylabel('y');
title('AB2 Solutions');

hold off

%% Explicit RK2
[x1,y1] = RK2(0,0.033,30);
[x2,y2] = RK2(0, 0.04,25);
[x3,y3] = RK2(0, 0.05,20);
[x4,y4] = RK2(0, 0.02,50);

figure
hold on
grid on

x = 0:0.0001:1;
exact = (50*2501^(1/2)*cos(- atan(1/50) + x))/2501 - (2500*exp(-50.*x))/2501;

plot(x3,y3, 'Color',[0 0.4470 0.7410]);
plot(x2,y2, 'Color', [0.8500 0.3250 0.0980]);
plot(x1,y1,'Color', [0.4660 0.6740 0.1880]);
plot(x4,y4,'Color', [0.9290 0.6940 0.1250]);

plot(x,exact, 'Color', [0.4940 0.1840 0.5560]);

legend('h = 0.05, 20 steps','h = 0.04, 25 steps', 'h = 0.033, 30 steps', 'h = 0.02, 50 steps', 'Exact Solution');
xlim([0 1]);
ylim([0 1.8]);
xlabel('x');
ylabel('y');
title('RK2 Solutions');

hold off

%% Explicit RK4

[x1,y1] = RK4(0,0.033,30);
[x2,y2] = RK4(0, 0.04,25);
[x3,y3] = RK4(0, 0.05,20);
[x4,y4] = RK4(0, 0.02,50);

figure
hold on
grid on

x = 0:0.0001:1;
exact = (50*2501^(1/2)*cos(- atan(1/50) + x))/2501 - (2500*exp(-50.*x))/2501;

plot(x3,y3, 'Color',[0 0.4470 0.7410]);
plot(x2,y2, 'Color', [0.8500 0.3250 0.0980]);
plot(x1,y1,'Color', [0.4660 0.6740 0.1880]);
plot(x4,y4,'Color', [0.9290 0.6940 0.1250]);

plot(x,exact, 'Color', [0.4940 0.1840 0.5560]);

legend('h = 0.05, 20 steps','h = 0.04, 25 steps', 'h = 0.033, 30 steps', 'h = 0.02, 50 steps', 'Exact Solution');
xlim([0 1]);
ylim([0 1.8]);
xlabel('x');
ylabel('y');
title('RK4 Solutions');

hold off


%% Global Error

exact = (50*2501^(1/2)*cos(- atan(1/50) + 1))/2501 - (2500*exp(-50.*1))/2501;

h = [2.^(-(6:10))];

feuler_e = zeros(1,length(h));
beuler_e = zeros(1,length(h));
trap_e = zeros(1,length(h));
midpoint_e = zeros(1,length(h));
AB2_e = zeros(1,length(h));
RK2_e = zeros(1,length(h));
RK4_e = zeros(1,length(h));

for i = 1:length(h)
    N = 1/h(i);
    y1 = 0;
    [x,y,Mf] = feuler(y1,h(i),N);
    feuler_e(i) = y(end)-exact;
    
    [x,y,Mb] = beuler(y1,h(i),N);
    beuler_e(i) = y(end)-exact;
    
    [x,y,Mt] = trapezoid(y1,h(i),N);
    trap_e(i) = y(end)-exact;
    
    [x,y,Mp] = midpoint(y1,h(i),N);
    midpoint_e(i) = y(end)-exact;
    
    [x,y,Ma] = AB2(y1,h(i),N);
    AB2_e(i) = y(end)-exact;
    
    [x,y,Mr2] = RK2(y1,h(i),N);
    RK2_e(i) = y(end)-exact;
    
    [x,y,Mr4] = RK4(y1,h(i),N);
    RK4_e(i) = y(end)-exact;
end



figure 

loglog(1./h, abs(feuler_e));
hold on
grid on
loglog(1./h, abs(beuler_e));
loglog(1./h, abs(trap_e));
loglog(1./h, abs(midpoint_e));
loglog(1./h, abs(AB2_e));
loglog(1./h, abs(RK2_e));
loglog(1./h, abs(RK4_e));


legend('Explicit Euler', 'Implicit Euler', 'Trapezoidal', 'Midpoint', 'AB2', 'RK2', 'RK4');
xlabel('1/h')
ylabel('|e| (absolute value of error)');
title('Global Error of Numerical Methods');
hold off

M1 = [Mf(1) Mb(1) Mt(1) Mp(1) Ma(1) Mr2(1) Mr4(1)]
e1 = [feuler_e(1) beuler_e(1) trap_e(1) midpoint_e(1) AB2_e(1) RK2_e(1) RK4_e(1)]

figure
plot(e,M);
grid on
hold on

title('Error vs Work');
xlabel('M work');
ylabel('error');

%% Problem 2: Eigenvalues

N = 10;
Tn = diag(-2*ones(1,N))+diag(1*ones(1,N-1),1)+diag(1*ones(1,N-1),-1);

Tneig = eig(Tn);
Tneig = sort(Tneig);

echeck = zeros(1,N);

for i =1:N
    
    values = -2*(1-cos(pi*i/(N+1)));
    echeck(i) = values;
    echeck = sort(echeck);
end
Tneig = Tneig.'
echeck = echeck
meig = zeros(1,20);

for N = 1:20
    
    Tp = diag(-2*ones(1,N))+diag(1*ones(1,N-1),1)+diag(1*ones(1,N-1),-1);
    eigp = eig(Tp);
    abseig = abs(eigp);
    meig(N) = max(abseig);
    
end

figure
hold on
grid on

N = 1:20;
checking = -2*(1-cos(pi.*N./(N+1)));
scatter(N, meig, 'filled');
plot(N, abs(checking));
title('Absolute Max Eigenvalue for Tn(NxN) vs. N');
xlabel('N');
ylabel('Max Eigen');
legend('Max Eigenvalue', 'lambda = abs(-2*(1-cos(pi*N/(N+1)))');


%% Functions

function[x,y, M] = feuler(y1,h,N)
    x = zeros(1,N);
    y = zeros(1,N);
    M = 0;
    y(1) = y1;
    for n = 1:N
        y(n+1) = y(n) + h*(-50*(y(n)-cos(x(n))));
        M=M+1;
        x(n+1) = n*h;   
    end
end

function[x,y,M] = beuler(y1,h,N)
    x = zeros(1,N);
    y = zeros(1,N);
    y(1) = y1;
    M=0;
    for n = 1:N
        x(n+1) = n*h;
        y(n+1) = (y(n) + h*(50*(cos(x(n+1)))))/(1+50*h);
        M=M+1;
        
    end        

end

  function [x,y,M] = midpoint(y1,h,N)
    x = zeros(1,N);
    y = zeros(1,N);
    y(1) = y1;
    M=0;
    for n = 1:N
        y(n+1) = y(n) + h*(-50*(y(n)-cos(x(n)+h/2)));
        M=M+1;
        x(n+1) = n*h;   
    end

  end
  
  function [x,y,M] = trapezoid(y1,h,N)
  
    x = zeros(1,N);
    y = zeros(1,N);
    y(1) = y1; 
    x(1) = 0;
    M=0;
    
    for n = 1:N
        x(n+1) = n*h;
        y(n+1) = (y(n) + (h)*(-25*y(n)+25*cos(x(n))+25*cos(x(n+1))))/(1+25*h);
        M=M+1;
    end  
    
  end
  
  function [x,y,M] = AB2(y1,h,N)
    x(1) = 0;
    y(:,1) = y1;
    M =0;
    f_ode = @(x,y) -50*(y-cos(x));
    k = 1;
    fValue =  f_ode( x(k), y(:,k) );
    xhalf = x(k) + 0.5 * h;
    yhalf = y(:,k) + 0.5 * h * fValue;
    fValuehalf = f_ode( xhalf, yhalf );

    x(1,k+1) = x(1,k) + h;
    y(:,k+1) = y(:,k) + h * fValuehalf;

    for k = 2 : N
        fValueold=fValue;
        fValue = f_ode( x(k), y(:,k) );
        M = M+1;
        x(1,k+1) = x(1,k) + h;
        y(:,k+1) = y(:,k) + h * ( 3 * fValue - fValueold ) / 2;
        M=M+1;
    end
    
  end
  
    function [x, y, M] = RK2(y1,h,N)
  
    x = zeros(1,N);
    y = zeros(1,N);
    y(1) = y1; 
    x(1) = 0;
    M=0;
    for n = 1:N
      k1 = -50*(y(n)-cos(x(n)));
      M=M+1;
      k2 = -50*(y(n)+h*k1/2 - cos(x(n)+h/2));
      M=M+1;
      y(n+1) = y(n) + h/6*(k1+2*k2); 
      x(n+1) = n*h;
    end  
  
  end
  
  function [x, y,M] = RK4(y1,h,N)
  
    x = zeros(1,N);
    y = zeros(1,N);
    y(1) = y1; 
    x(1) = 0;
    M = 0;
    for n = 1:N
      k1 = -50*(y(n)-cos(x(n)));
      M=M+1;
      k2 = -50*(y(n)+h*k1/2 - cos(x(n)+h/2));
      M=M+1;
      k3 = -50*(y(n)+h*k2/2 - cos(x(n)+h/2));
      M=M+1;
      k4 = -50*(y(n)+h*k3 - cos(x(n)+h));
      M=M+1;
      y(n+1) = y(n) + h*(1/6*k1+2/6*k2+2/6*k3+1/6*k4); 
      x(n+1) = n*h;
    end  
  
  end
  