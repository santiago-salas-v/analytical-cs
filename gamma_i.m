function gamma_i
theta_j     = reshape(theta_j,1,numel(theta_j));
s_j         = reshape(s_j,1,numel(s_j));
gamma_i_C   = exp(1 - J_i + log(J_i)-5*q_i*(1-J_i/L_i+log(J_i/L_i)));
gamma_i_R   = exp(q_i*(1-log(s_i)-sum(theta_j./s_j*tau_ij)));
end