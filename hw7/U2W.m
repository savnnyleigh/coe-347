function [W] = U2W(U,g)
W = U*0;

W(1) = U(1);
W(2) = U(2)/U(1);
W(3) = (g-1)*(U(3)-U(2)^2/(2*U(1)));

end