function kl = KL_Univar_symmetric(p_mean, p_var, q_mean, q_var)

% calculate KL divergence for two univariate gaussian probability
% distributions

% A) p vs q
kl_pq = [];

a = log(sqrt(q_var)/sqrt(p_var));
b = p_var + power((p_mean - q_mean),2);
c = 2 * q_var;

kl_pq = a + (b./c) - 0.5;

% B) q vs p
kl_qp = [];

a = log(sqrt(p_var)/sqrt(q_var));
b = q_var + power((q_mean - p_mean),2);
c = 2 * p_var;

kl_qp = a + (b./c) - 0.5;

% kl symmetric = result from A vs result from B
kl = (kl_pq + kl_qp)/2;