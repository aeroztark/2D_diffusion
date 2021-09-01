
function[A] = ConstCoeffImplicitCNMatrix(k,dxz,i,j)
% This function generates the matrix for solution of 2D diffusion equation
% by Direct implicit method using Crank-Nicholson scheme.
% This function works for uniform grid (dx=dz) and constant diffusion
% coefficient.

% k -> diffusion constant
% dxz -> grid size for uniform grid (dx = dz)
% i -> number of x values in u(length of x dimension)
% j -> number of z values in u (length of z dimension)

a = 1+(2*k/(dxz^2));
c = -k/(2*dxz^2);
ij = i*j;

main_diag = a.*ones(ij,1);
plusi_diag = c.*ones(ij,1);
mini_diag = plusi_diag;

c_unit = c.*ones(i,1);
min1c_unit = c_unit;
min1c_unit(end) = 0;
min1_diag = repmat(min1c_unit,j,1);

plus1c_unit = c_unit;
plus1c_unit(1) = 0;
plus1_diag = repmat(plus1c_unit,j,1);

A = spdiags([mini_diag,min1_diag,main_diag,plus1_diag,plusi_diag],[-i,-1,0,1,i],ij,ij);

end
