function[u_new] = ImplicitMethod(k,u_old,dxy,i,j)
% Thus function solves the 2D diffusion equation using a Direct Implicit Method based on Cranck-Nocholson scheme.

% created by: Sarthak Srivastava, 01-Sep-2021

        % get coefficient matrix
        A = ConstCoeffImplicitCNMatrix(k,dxy,i,j);
        % solve the linear system
        b = reshape(u_old(:,:)',[],1);
        u_new = A\b; % Direct method using LU factorization
        % reshape the solution as matrix
        u_new = reshape(u_new,[i,j])';
        % apply BCs
        u_new(1,:) = 0;
        u_new(end,:) = 0;
        u_new(:,1) = 0;
        u_new(:,end) = 0;    
end