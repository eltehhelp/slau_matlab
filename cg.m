function [x,k,c] = cg(A, b, x, tol, maxit)
% CG Solves Ax=b using Conjugate Gradient
%
% Inputs:
%     A     Matrix to solve in systemm
%     b     Right hand side vector
%     x     Initial guess for iterative method
%     tol   Desired tolerance to achieve in the solution
%     maxit Maximum number of iterations
%
% Outputs:
%     x  Solution
%     k  Number of iterations taken
%     c  Vector of residual ||b-A*x_k|| at each step

nb = norm(b);

c = zeros(maxit,1);

r = b - A * x;
p = r;

rr = r'*r;
k = 1;
c(k) = sqrt(rr);

while (c(k)/nb > tol) && (k < maxit)
    Ap     = A*p;
    alpha  = rr/(p'*Ap);
    x      = x + alpha*p;
    r      = r - alpha*Ap;
    
    RR     = rr;
    rr     = r'*r;
    beta   = rr/RR;
    p      = r + beta*p;

    k = k + 1;
    c(k) = sqrt(rr);
end
c = c(1:k);

end
