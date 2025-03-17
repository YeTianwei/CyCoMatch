function C_opt = admm_optimization(C0, beta, mu, G, cycles, rho, max_iterations, tolerance)

    n = size(C0{1,1}, 1);
    C = C0;  
    Z = C0;  
    
    U = cell(size(C0));
    for i = 1:size(C0)
        for j = 1:size(C0)
            if ~isempty(C{i,j})
                U{i, j} = zeros(size(C{i, j}));
            end
        end
    end
    
    [~, edges] = adj_mat_2_edge(G);
    edges = edges';
    
    epsilon = 1e-6;
    
    for iteration = 1:max_iterations
               
        C_old = C;
        Z_old = Z;

        C = update_C(C_old, Z, U, rho, C0, edges);

        Z = update_Z(Z_old, C, U, beta, mu, edges, cycles, n, epsilon);

        U = update_U(U, C, Z, edges);

        if check_convergence(C, C_old, Z, Z_old, tolerance)
            break;
        end
    end

    C_opt = C;
end

function C_new = update_C(C_old, Z, U, rho, C0, edges)

    C_new = C_old;
    for k = 1:size(edges, 1)
        i = edges(k, 1);
        j = edges(k, 2);
        C_new{i, j} = (C0{i, j} + rho * (Z{i, j} - U{i, j})) / (1 + rho);
        
        C_new{j, i} = (C0{j, i} + rho * (Z{j, i} - U{j, i})) / (1 + rho);
        
    end
    
end

function Z_new = update_Z(Z_old, C, U, beta, mu, edges, cycles, n, epsilon)

    Z_new = Z_old;
    Z1 = cell(length(Z_old));
    Z2 = cell(length(Z_old));
    Z3 = cell(length(Z_old));
    
    for k = 1:size(edges, 1)
        i = edges(k, 1);
        j = edges(k, 2);
        
        Z1{i,j} = C{i,j} + U{i,j};
        Z1{j,i} = C{j,i} + U{j,i};
        
    end

    sum_weight = 0;

    % cycle constriant
    for c = 1:length(cycles)
        cycle = cycles{c}.cycle;
        
        for k = 1:length(cycle)
            prod = eye(n);
            i = cycle(k);
            
            if k==2
                j = cycle(k+1);
            else
                j = cycle(mod(k+1,length(cycle))); 
            end
            
            if k == 1
                z = cycle(k+2); 
            else
                z = cycle(mod(k+2,length(cycle)));
            end
            
            prod = prod * Z_old{i, j} * Z_old{j,z};
            Z2{i,j} = beta * cycles{c}.weight * (prod - eye(size(prod)));           
        end
        sum_weight = sum_weight + cycles{c}.weight;
    end
        
    for k = 1:size(edges, 1)
        i = edges(k, 1);
        j = edges(k, 2);
        
        Z1{i,j} = C{i,j} + U{i,j};
        Z1{j,i} = C{j,i} + U{j,i};
        
        Z3{i,j} = mu * Z_old{i,j} * (Z_old{i,j}' * Z_old{i,j} - eye(size(Z_old{i,j}))); 
        Z3{j,i} = mu * Z_old{j,i} * (Z_old{j,i}' * Z_old{j,i} - eye(size(Z_old{j,i}))); 
        
        Z3{i,j} = mu * Z_old{i,j} * (Z_old{i,j}' * Z_old{i,j} - eye(size(Z_old{i,j}))) + epsilon * eye(size(Z_old{i,j}));
        
        Z_new{i,j} = (Z1{i,j} + Z2{i,j} + Z3{i,j}) / (sum_weight  + beta + mu);
        Z_new{j,i} = (Z1{j,i} + Z2{j,i} + Z3{j,i}) / (sum_weight  + beta + mu);

    end
       
end

function U_new = update_U(U, C, Z, edges)
    
    U_new = U;
    for k = 1:size(edges, 1)
        i = edges(k, 1);
        j = edges(k, 2);
        U_new{i, j} = U{i, j} + C{i, j} - Z{i, j};
        U_new{j, i} = U{j, i} + C{j, i} - Z{j, i};
    end
end

function converged = check_convergence(C, C_old, Z, Z_old, tolerance)

    diff_C = 0;
    diff_Z = 0;
    for i = 1:size(C)
        for j = 1:size(C)
            diff_C = diff_C + norm(C{i,j} - C_old{i,j}, 'fro');
            diff_Z = diff_Z + norm(Z{i,j} - Z_old{i,j}, 'fro');
        end
    end
    converged = (diff_C < tolerance) && (diff_Z < tolerance);
end


% Compute edge matrix
function [edgeMatrix, edges] = adj_mat_2_edge(A)

    numV = size(A, 1);
    [rows, cols, ~] = find(A);
    ids = find(rows < cols);
    rows = rows(ids);
    cols = cols(ids);

    edges = [rows, cols]';
    numE = size(edges, 2);
    edgeMatrix = zeros(numV, numV);
    for eId = 1 : numE
        edgeMatrix(edges(1, eId), edges(2, eId)) = eId;
        edgeMatrix(edges(2, eId), edges(1, eId)) = eId;
    end
    
end
