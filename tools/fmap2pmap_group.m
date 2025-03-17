
function T = fmap2pmap_group(shapes,C)

T = cell(length(C));

for i = 1:length(C)
    for j = 1:length(C)
        
        T{i,j} = knnsearch(gpuArray(shapes{j}.evecs * C{i,j}'),gpuArray(shapes{i}.evecs)); % gpu version
        T{i,j} = gather(T{i,j});
       
        % T{i,j} = knnsearch(shapes{j}.evecs * C{i,j}',shapes{i}.evecs); % cpu version

        T{i,j} = fast_pMap_NNinterp(T{i,j},shapes{i});
        T{i,j} = T{i,j}(:);

    end
end