T=2*2*pi; 
exact = load('exact.dat');
n = [25, 50, 100, 200, 400];
h = 2*pi./(n+1);
% Calculating FOU errors

% sigma = 0.25

efou25 = zeros(1,5);
r1 = load('forward-upwind_n25_sigma0.250000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r1(i,end)-sin(r1(i,1))).^2;
    y = y + r;
end
efou25(1) = 1/25*sqrt(y);

r2 = load('forward-upwind_n50_sigma0.250000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r2(i,end)-sin(r2(i,1))).^2;
    y = y + r;
end
efou25(2) = 1/50*sqrt(y);

r3 = load('forward-upwind_n100_sigma0.250000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r3(i,end)-sin(r3(i,1))).^2;
    y = y + r;
end
efou25(3) = 1/100*sqrt(y);

r4 = load('forward-upwind_n200_sigma0.250000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r4(i,end)-sin(r4(i,1))).^2;
    y = y + r;
end
efou25(4) = 1/200*sqrt(y);

r5 = load('forward-upwind_n400_sigma0.250000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r5(i,end)-sin(r5(i,1))).^2;
    y = y + r;
end
efou25(5) = 1/400*sqrt(y);

% sigma = 0.5

efou5 = zeros(1,5);
r1 = load('forward-upwind_n25_sigma0.500000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r1(i,end)-sin(r1(i,1))).^2;
    y = y + r;
end
efou5(1) = 1/25*sqrt(y);

r2 = load('forward-upwind_n50_sigma0.500000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r2(i,end)-sin(r2(i,1))).^2;
    y = y + r;
end
efou5(2) = 1/50*sqrt(y);

r3 = load('forward-upwind_n100_sigma0.500000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r3(i,end)-sin(r3(i,1))).^2;
    y = y + r;
end
efou5(3) = 1/100*sqrt(y);

r4 = load('forward-upwind_n200_sigma0.500000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r4(i,end)-sin(r4(i,1))).^2;
    y = y + r;
end
efou5(4) = 1/200*sqrt(y);

r5 = load('forward-upwind_n400_sigma0.500000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r5(i,end)-sin(r5(i,1))).^2;
    y = y + r;
end
efou5(5) = 1/400*sqrt(y);


% sigma = 0.75

efou75 = zeros(1,5);
r1 = load('forward-upwind_n25_sigma0.750000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r1(i,end)-sin(r1(i,1))).^2;
    y = y + r;
end
efou75(1) = 1/25*sqrt(y);

r2 = load('forward-upwind_n50_sigma0.750000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r2(i,end)-sin(r2(i,1))).^2;
    y = y + r;
end
efou75(2) = 1/50*sqrt(y);

r3 = load('forward-upwind_n100_sigma0.750000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r3(i,end)-sin(r3(i,1))).^2;
    y = y + r;
end
efou75(3) = 1/100*sqrt(y);

r4 = load('forward-upwind_n200_sigma0.750000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r4(i,end)-sin(r4(i,1))).^2;
    y = y + r;
end
efou75(4) = 1/200*sqrt(y);

r5 = load('forward-upwind_n400_sigma0.750000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r5(i,end)-sin(r5(i,1))).^2;
    y = y + r;
end
efou75(5) = 1/400*sqrt(y);

% Calculating LW errors

% sigma = 0.25

elw25 = zeros(1,5);
r1 = load('lax-wendroff_n25_sigma0.250000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r1(i,end)-sin(r1(i,1))).^2;
    y = y + r;
end
elw25(1) = 1/25*sqrt(y);

r2 = load('lax-wendroff_n50_sigma0.250000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r2(i,end)-sin(r2(i,1))).^2;
    y = y + r;
end
elw25(2) = 1/50*sqrt(y);

r3 = load('lax-wendroff_n100_sigma0.250000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r3(i,end)-sin(r3(i,1))).^2;
    y = y + r;
end
elw25(3) = 1/100*sqrt(y);

r4 = load('lax-wendroff_n200_sigma0.250000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r4(i,end)-sin(r4(i,1))).^2;
    y = y + r;
end
elw25(4) = 1/200*sqrt(y);

r5 = load('lax-wendroff_n400_sigma0.250000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r5(i,end)-sin(r5(i,1))).^2;
    y = y + r;
end
elw25(5) = 1/400*sqrt(y);

% sigma = 0.5

elw5 = zeros(1,5);
r1 = load('lax-wendroff_n25_sigma0.500000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r1(i,end)-sin(r1(i,1))).^2;
    y = y + r;
end
elw5(1) = 1/25*sqrt(y);

r2 = load('lax-wendroff_n50_sigma0.500000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r2(i,end)-sin(r2(i,1))).^2;
    y = y + r;
end
elw5(2) = 1/50*sqrt(y);

r3 = load('lax-wendroff_n100_sigma0.500000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r3(i,end)-sin(r3(i,1))).^2;
    y = y + r;
end
elw5(3) = 1/100*sqrt(y);

r4 = load('lax-wendroff_n200_sigma0.500000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r4(i,end)-sin(r4(i,1))).^2;
    y = y + r;
end
elw5(4) = 1/200*sqrt(y);

r5 = load('lax-wendroff_n400_sigma0.500000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r5(i,end)-sin(r5(i,1))).^2;
    y = y + r;
end
elw5(5) = 1/400*sqrt(y);


% sigma = 0.75

elw75 = zeros(1,5);
r1 = load('lax-wendroff_n25_sigma0.750000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r1(i,end)-sin(r1(i,1))).^2;
    y = y + r;
end
elw75(1) = 1/25*sqrt(y);

r2 = load('lax-wendroff_n50_sigma0.750000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r2(i,end)-sin(r2(i,1))).^2;
    y = y + r;
end
elw75(2) = 1/50*sqrt(y);

r3 = load('lax-wendroff_n100_sigma0.750000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r3(i,end)-sin(r3(i,1))).^2;
    y = y + r;
end
elw75(3) = 1/100*sqrt(y);

r4 = load('lax-wendroff_n200_sigma0.750000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r4(i,end)-sin(r4(i,1))).^2;
    y = y + r;
end
elw75(4) = 1/200*sqrt(y);

r5 = load('lax-wendroff_n400_sigma0.750000.dat');
y = 0; 
for i = 1:length(r1)
    r = (r5(i,end)-sin(r5(i,1))).^2;
    y = y + r;
end
elw75(5) = 1/400*sqrt(y);

% Plot errors

figure


loglog(1./h, efou25);
grid on
hold on
loglog(1./h, efou5);
loglog(1./h, efou75);
loglog(1./h, elw25, '--');
loglog(1./h, elw5, '--');
loglog(1./h, elw75, '--');

xlabel('1/h');
ylabel('error');
title('Error vs 1/h');
legend('FOU sigma =0.25', 'FOU sigma =0.5', 'FOU sigma =0.75', 'LW sigma =0.25', 'LW sigma =0.5', 'LW sigma =0.75');