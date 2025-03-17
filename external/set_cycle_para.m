function Para = set_cycle_para()

% setting admm param to generate cycle consistency basis

Para = struct;

Para.mu_init = 1;
Para.mu_rho = 1.0060;

Para.exact_solver = 1;

Para.cg_eps = 1.0000e-06;
Para.cg_iters = 50;

Para.using_cvx = 0;

Para.num_iters = 600;

end