function  [T0, dim, C1] = LMD_Filter(shapes, Dist, T0, threshold)
    
    num_shapes = length(shapes);
    
    dim = zeros(num_shapes);
    
    C1 = cell(num_shapes);
    
    for i = 1:num_shapes
        
        for j = 1:num_shapes
            
            num_vert = min(shapes{i}.nv,shapes{j}.nv);
            good_vert = 1:num_vert;
            distortion = zeros(num_vert, 1);
            dist_source_matrix = Dist{i}(good_vert, good_vert);
            
            mesh_areas = full(diag(shapes{j}.A));
            dist_target_matrix = Dist{j}(:, T0{i, j}(good_vert));
            
            R_max = max(max(Dist{j}));
            
            for k = 1:num_vert
                idx_source = find(dist_source_matrix(:, k) ~= 0);
                dist_source = dist_source_matrix(idx_source, k);
                dist_target = dist_target_matrix(idx_source, k);
                dist_target(dist_target == 0) = R_max;
                
                max_distance = max(dist_source);
                
                distortion(k) = sum(((abs(dist_source - dist_target)) / max_distance) .* mesh_areas(idx_source)) / sum(mesh_areas(idx_source));
            end
            
            sub_landmarks = find(distortion < threshold);
            non_landmarks = setdiff(1:num_vert, sub_landmarks)';
            non_landmarks_corr = setdiff(1:num_vert, T0{i, j}(sub_landmarks))';
            
            H = shapes{i}.evecs(sub_landmarks, :)' * shapes{j}.evecs(T0{i, j}(sub_landmarks), :);
            [U1, D_matrices, V1] = svd(H);
            dim(i, j) = findK(diag(D_matrices));
            C1{i,j} = U1 * V1';
            
            projC = shapes{i}.evecs(sub_landmarks, 1:dim(i, j))' * shapes{j}.evecs(T0{i, j}(sub_landmarks), 1:dim(i, j));
            [U, ~, V] = svd(projC);
            C_fmap = U * V';
            
            spec1_sub = shapes{i}.evecs(non_landmarks, 1:dim(i, j));
            spec2_sub = shapes{j}.evecs(non_landmarks_corr, 1:dim(i, j));
            
            T_temp = knnsearch(gpuArray(spec2_sub * C_fmap'), gpuArray(spec1_sub)); % gpu version
            T_temp = gather(T_temp);

            % T_temp = knnsearch(spec2_sub * C_fmap', spec1_sub); % cpu version

            T_temp = T_temp(:);
            
            T0{i, j}(non_landmarks) = non_landmarks_corr(T_temp);
            
            T0{i,j} = fast_pMap_NNinterp(T0{i,j},shapes{i});
        end
    end

end