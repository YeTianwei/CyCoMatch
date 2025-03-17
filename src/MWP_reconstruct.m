function Cout = MWP_reconstruct(shapes, T0, Nf, num_iters, dim)

n = length(shapes);
Cout = cell(n);

for i = 1:n
    
    for j = 1:n
        T = T0{i,j};
        d = dim(i,j);
        
        g = gsp_design_meyer(shapes{i}.evals(d),Nf);
        ref_g = gsp_design_meyer(shapes{j}.evals(d),Nf);
        
        for it = 1:num_iters
            C = zeros(d);
            C_fmap = shapes{i}.evecs(:,1:d)\shapes{j}.evecs(T,1:d);

            for s=1:Nf
                ref_fs = sparse(1:d,1:d,ref_g{s}(shapes{j}.evals(1:d,:)));
                fs = sparse(1:d,1:d,g{s}(shapes{i}.evals(1:d,:)));
                C = C + fs*C_fmap*ref_fs;
            end
            
            T = knnsearch(gpuArray(shapes{j}.evecs(:,1:d)*C'),gpuArray(shapes{i}.evecs(:,1:d))); % gpu version
            T = gather(T);

            % T = knnsearch(shapes{j}.evecs(:,1:d)*C',shapes{i}.evecs(:,1:d)); % cpu version

            T = T(:);
            T = fast_pMap_NNinterp(T,shapes{i});
            
        end
        Cout{i,j} = C;
  
    end
end

end

