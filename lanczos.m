function [U, H] = lanczos(funMult, x, tol, maxIter)
% Perform m steps of the Lanczos process, storing the full basis U
% with reorthogonalization on Hermitian A with starting vector x.
%
% This code will malfunction if A is not square, if x is of different
% dimension than A, if x = 0, or if m is greater than the dimemion
% of A.

x = x / norm(x);

alpha_arr = [];
beta_arr = [];
n = length(x);
% we store all vectors
U = x;
beta = 0;
maxIter = min(n,maxIter);
for j = 1:maxIter
    z = funMult(U(:,end));
    alpha = U(:,end)'*z;
    u = z - U(:,end)*alpha;
    if j > 1
        u = u - U(:,end-1)*beta;
    end
    % Reorthognalization
    alphas = U'*u;
    u = u - U*alphas;
    alpha = alpha + alphas(end);
    alpha_arr(end+1) = alpha;
    
    beta = norm(u);
    if beta < tol
        break;
    end
    beta_arr(end+1) = beta;
    u = u / beta;
    U = [U,u];
end

H = diag(alpha_arr) + diag(beta_arr, 1) + diag(beta_arr, -1);

end