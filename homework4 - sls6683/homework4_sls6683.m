%% Problem 2

T = readtable('solutionA_N10000.dat');

xexact = T{:,1};
yexact = T{:,2};

[x10,U10] = central(10);

figure
hold on
grid on

plot(xexact,yexact);
plot(x10,U10, '-s');


xlim([0, 1]);

xlabel('x');
ylabel('u(x)');
title('u(x) estimations');
legend('Exact Solution', 'N=10');

hold off

%{
M = table(x10',U10');
writetable(M, 'Problem2data.xlsx');
%}
%% Problem 3

N = [5, 10, 20, 40, 80, 160, 320, 640, 1280];
h = [1./(N+1)];

E = zeros(9,1);
e = zeros(9,1);

l = @(x) -1.*cos((4*pi).*x)+ 2.*x +1;

for i = 1:9
    
    [xi,ui] = central(N(i));
    u = l(xi);
    for j = 1:N(i)+2
        p(j) = (ui(j)-u(j)).^2;
    end
    
    E(i) = (sum(p)).^(1/2);
    e(i) = 1./N(i).*(sum(p)).^(1/2);
end


figure

loglog(1./h, E);
hold on
grid on
loglog(1./h, e);

title('log-log plot of E and e vs 1/h');
xlabel('1/h');
ylabel('Error');
legend('E', 'e');

hold off

%% Problem 5
T = readtable('solutionB_N10000.dat');

xexact = T{:,1};
yexact = T{:,2};


[x10,u10, A, g,u, rhs] = oneside(10);


figure
hold on
grid on

plot(xexact,yexact, 'LineWidth', 2);
plot(x10,u10, '-s');

xlim([0, 1]);

xlabel('x');
ylabel('u(x)');
title('u(x) estimations');
legend('Exact Solution', 'N=10');

hold off

M = table(x10',u10');
writetable(M, 'Problem5data.xlsx');

%% Problem 6

N = [5, 10, 20, 40, 80, 160, 320, 640, 1280];
h = [1./(N+1)];

E = zeros(9,1);
e = zeros(9,1);

l = @(x) -1.*cos((4*pi).*x)+ 2.*x +1;

for i = 1:9
    
    [xi,ui] = oneside(N(i));
    u = l(xi);
    for j = 1:N(i)+2
        p(j) = (ui(j)-u(j)).^2;
    end
    
    E(i) = (sum(p)).^(1/2);
    e(i) = 1./N(i).*(sum(p)).^(1/2);
end


figure

loglog(1./h, E);
hold on
grid on
loglog(1./h, e);

title('log-log plot of E and e vs 1/h');
xlabel('1/h');
ylabel('Error');
legend('E', 'e');

hold off


%% Functions

function [xi, Ui] = central(n);
    N = n;
    h = 1./(N+1);
    xi = 0:h:1;

    f = @(x) (4*pi)^2.*cos((4*pi).*x);
    t = f(xi)';

    A = diag((-2*ones(1,N)))+diag((1*ones(1,N-1)),1)+diag((1*ones(1,N-1)),-1);
    u = zeros(N,1);
    g = zeros(N,1);

    a = 0;
    b = 2;
    g(1) = a;
    g(N) = -b;

    rhs = t(2:N+1).*h.^2+g;
    u = A\rhs;

    Ui = [0 u' 2];

end

function [xi, Ui, A, g, u, rhs] = oneside(n);
    N = n+1;
    h = 1./(N);
    xi = 0:h:1;

    f = @(x) (4*pi)^2.*cos((4*pi).*x);
    t = f(xi)';
    t(1) = 0;

    A = diag((-2*ones(1,N)))+diag((1*ones(1,N-1)),1)+diag((1*ones(1,N-1)),-1);
    A(1,:) = [(-3/2) 2 (-1/2) zeros(1,N-3)];
    u = zeros(N,1);
    g = zeros(N,1);

    b = 2;
    
    g(1) = 10.*h;
    g(N) = -b;

    rhs = t(1:N).*h.^2+g;
    u = A\rhs;

    Ui = [u' 2];

end