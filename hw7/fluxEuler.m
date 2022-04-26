function [F] = fluxEuler(U,g)

W = U2W(U,g);

%F(1) = W(2)*U(1); same as below
F(1) = U(2);

F(2) = W(2)*U(2)+W(3);

F(3) = W(2)*(U(3)+W(3));

end