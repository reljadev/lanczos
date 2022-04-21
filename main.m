% ----------------------- TIME AND PRECISION ----------------------------

sizes = [50, 100, 150, 200, 250, 300];
eps = 1e-8;
global Ua Ub Sa Sb Va Vb;
    
for ii = 1:length(sizes)
    
    % create matrices A and B
    M = sizes(ii);
    N = sizes(ii);
    x = 0.1:0.1:M;
    y = 0.1:0.1:N;
    [X, Y] = meshgrid(x, y);
    A = 1 ./ (X + Y);
    B = 1 ./ sqrt(X.^2 + Y.^2);
    C = A .* B;
    
    % calculate svd of A and B
    [Ua, Sa, Va] = truncated_svd(A, eps);
    [Ub, Sb, Vb] = truncated_svd(B, eps);
    
    % calculate svd of C directly
    t = tic;
    [Uc, Sc, Vc] = truncated_svd(C, eps);
    time_svd(ii) = toc(t);
    
    svd_errors(ii) = norm(C - Uc * Sc * Vc');
    
    % calculate approximation of C by Lanczos + fast mv
    % initialize random vector
    x = rand(10*N,1);
    % || C*C' - Q*H*Q' || < tol
    t = tic;
    [Q, H] = lanczos(@mvMult_times, x, 1e-8, 30);
    
    % || H - V*L*V' || < tol
    [V, L] = eig(H);
    
    % || C*C' - Uc*L*Uc' || < 2*tol
    Uc = Q * V;
    % get sigmas, s_i = sqrt(l_i)
    sigmas = sqrt(diag(L));
    % v_i = C'*u_i / s_i
    Vc = [];
    m = size(Uc, 2);
    for j = 1:m
        u = Uc(:, j);
        s = sigmas(j);
        %v = (C' * u) ./ s;
        v = (mvMult_transpose(u)) ./ s;
        Vc = [Vc, v];
    end
    % Sc = sqrt(L)
    Sc = diag(sigmas);
    
    time_approx(ii) = toc(t);
    
    % error of our approximation
    approx_errors(ii) = norm(C - Uc * Sc * Vc', 'fro');
    
end

% plot times
figure()
plot(sizes, time_svd, 'bx-', sizes, time_approx, 'gx-')
set(gca,'fontsize',10)
xlabel('matrix size')
ylabel('time')
legend({'truncated SVD', 'lanczos + fast mvMult'}, 'Location', 'NorthWest');
grid on;

% plot error of approximations
figure()
semilogy(sizes, svd_errors, 'bx-', sizes, approx_errors, 'gx-')
set(gca,'fontsize',10)
xlabel('matrix size')
ylabel('error')
legend({'truncated SVD', 'lanczos + fast mvMult'}, 'Location', 'NorthWest');
grid on;

% ------------------------ RANK COMPARISON -----------------------------

eps = 1e-4;
    
% create matrices A and B
M = 1;
N = 2;
x = 0.1:0.1:M;
y = 0.1:0.1:N;
[X, Y] = meshgrid(x, y);
A = 1 ./ (X + Y);
B = 1 ./ sqrt(X.^2 + Y.^2);
C = A .* B;
    
% calculate svd of A and B
[Ua, Sa, Va] = truncated_svd(A, eps);
[Ub, Sb, Vb] = truncated_svd(B, eps);
    
% calculate svd of C directly
[Uc, Sc, Vc] = truncated_svd(C, eps);
    
% rank and accuracy of SVD and Hadamard representation
rankC = rank(Uc * Sc * Vc');
errorC = norm(C - Uc * Sc * Vc','fro');
rankC2 = rank(kr(Ua', Ub')' * kron(Sa, Sb) * kr(Va', Vb'));
errorC2 = norm(C - kr(Ua', Ub')' * kron(Sa, Sb) * kr(Va', Vb'),'fro');
    
    
% calculate approximation of C by Lanczos + fast mv
% initialize random vector
x = rand(10*N,1);
% || C*C' - Q*H*Q' || < tol
[Q, H] = lanczos(@mvMult_times, x, 1e-8, 30);
    
% || H - V*L*V' || < tol
[V, L] = eig(H);
    
% || C*C' - Uc*L*Uc' || < 2*tol
Uc = Q * V;
% get sigmas, s_i = sqrt(l_i)
sigmas = sqrt(diag(L));
% v_i = C'*u_i / s_i
Vc = [];
m = size(Uc, 2);
for j = 1:m
    u = Uc(:, j);
    s = sigmas(j);
    %v = (C' * u) ./ s;
    v = (mvMult_transpose(u)) ./ s;
    Vc = [Vc, v];
end
% Sc = sqrt(L)
Sc = diag(sigmas);
    
time_approx(ii) = toc(t);
    
% error of our approximation
approx_errors(ii) = norm(C - Uc * Sc * Vc', 'fro');
    
% calculate rank and accuracy of our approximation
rankC3 = rank(Uc * Sc * Vc');
errorC3 = norm(C - Uc * Sc * Vc','fro');

fprintf("C truncSVD rank: %d\n", rankC);
fprintf("C approximation rank: %d\n", rankC2);
fprintf("C our approx rank: %d\n", rankC3)
fprintf("Error of truncSVD: %f\n", errorC)
fprintf("Error of approximation: %f\n", errorC2)
fprintf("Error of our approximation: %f\n", errorC3)
