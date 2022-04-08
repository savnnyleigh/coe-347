

syms h

%FOU |g| with h = h == LW |g| with h = 0.1
equality =  (1-.5+.5*cos(h))^2+(.5*sin(h))^2 == (1-.5^2+.5^2*cos(.1))^2+(.5*sin(.1))^2; 
double(solve(equality,h))