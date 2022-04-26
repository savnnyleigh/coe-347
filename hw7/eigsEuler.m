function [lambda] = eigsEuler(U,g)

W = U2W(U,g);

a = sqrt(g*W(3)/W(1));

lambda = [W(2)-a W(2) W(2)+a];

end