function kl = KL_Multivar_symmetric(p_mean, p_cov, q_mean, q_cov, epsilon)

% calculate KL divergence for two multivariate gaussian probability
% distributions
% mean is a vector
% covar is the full covariance matrix

% epsilon used for ridge regression ('Tikhonov Regularization') to deal
% with singular numerically unstable covariance matrices
% A(-1) * x |-becomes-> [[A(T) * A] + epsilon^2 * eye(sizeA)] \ A(T) * x

assert(max(size(p_mean))==max(size(q_mean)),'p and q have different dimensions');

%%%%%%%%% should not need ridge regression if the shrinkage estimation
%%%%%%%%% works appropriately
% A) p vs q
kl_pq = [];

a = log(det(q_cov)/det(p_cov));
b = inv(q_cov) * p_cov;
c = q_mean - p_mean;
d = inv(q_cov) * c;
N = size(p_mean,1);

kl_pq = 0.5 * [a + trace(b) + (c' * d) - N];

% B) q vs p
kl_qp = [];

a = log(det(p_cov)/det(q_cov));
b = inv(p_cov) * q_cov;
c = p_mean - q_mean;
d = inv(p_cov) * c;
N = size(q_mean,1);

kl_qp = 0.5 * [a + trace(b) + (c' * d) - N];

% kl symmetric = result from A vs result from B
kl = (kl_pq + kl_qp)/2;

%%%%%%%%%%%% using ridge regression
%%%% still had issue of unstable covariance estimate to do determinant...
%%%% and overall should not be using a poorly conditioned covariance matrix
%%%% even if i could make it work

% % A) p vs q
% kl_pq = [];
% 
% a = log(det(q_cov)/det(p_cov));
% b = RidgeRegression(q_cov,p_cov,epsilon);
% c = q_mean - p_mean;
% d = RidgeRegression(q_cov,c,epsilon);
% N = size(p_mean);
% 
% kl_pq = 0.5 * [a + trace(b) + (c' * d) - N];
% 
% % B) q vs p
% kl_qp = [];
% 
% a = log(det(p_cov)/det(q_cov));
% b = RidgeRegression(p_cov,q_cov,epsilon);
% c = p_mean - q_mean;
% d = RidgeRegression(p_cov,c,epsilon);
% N = size(q_mean);
% 
% kl_qp = 0.5 * [a + trace(b) + (c' * d) - N];
% 
% % kl symmetric = result from A vs result from B
% kl = (kl_pq + kl_qp)/2;