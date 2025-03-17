function [C_opt,T] = Multi_CCB(shapes, T0, D, G)

    nshapes = length(shapes); T = cell(nshapes);
    
    % MWP parameter
    Nf = 6; num_iters = 7;  
    
    % add functional maps parameter
    alpha = 0.3;
    
    % ADMM parameter
    beta = 1; mu = 1; rho = 1e-3; max_iter = 100; tol = 1e-4;
    
    % compute cycles numbers
    num_cycles = (nshapes * (nshapes - 1) * (nshapes - 2)) / 2;
    
    % set parameter to generate cycles
    Para_cycle = set_cycle_para();
    
    % generate cycles from graph G
    cycles = cycle_basis_generator(G, num_cycles, Para_cycle);
    
    for it = 1:2
        
        fprintf('Iteration %d\n', it)
        
        threshold = [0.12,0.1];
        
        fprintf('Identify landmarks...\n')
        
        [TMap,dim,C1] = LMD_Filter(shapes, D, T0, threshold(it));
        
        fprintf('Reconstruct functional maps...\n')
        
        C2 = MWP_reconstruct(shapes, TMap, Nf, num_iters, dim);
        
        C = add_fmap(C1, C2, alpha);
        
        fprintf('Cycle Consistency Refinement...\n')
        
        C_opt = admm_optimization(C, beta, mu, G, cycles, rho, max_iter, tol);
    
        T = fmap2pmap_group(shapes, C_opt);
        
        T0 = T;
        
    end

end