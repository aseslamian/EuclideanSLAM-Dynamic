% function S = covariance_root( C )
% 
% EPS = 10^(-5);
% ZERO = 10^(-8);
% [~,error] = cholcov(C,0);
% 
% if error == 0
%     S = (chol(C))'; 
% else
%     [V, D] = eig(C); % eigen value decomposition
%     D(D <= ZERO) = EPS; % correcting eigen values if be less than zero
%     S = (chol(V*D*V'))';
% end
%     
% end


function S = covariance_root(C)

% Singular Value Decomposition
[U, Sigma, V] = svd(C);

% Compute the square root of the diagonal matrix
Sigma_root = diag(sqrt(diag(Sigma)));

% Calculate the matrix square root
S = U * Sigma_root * V';

end
