function D = calc_dist_mat(vert,tri)

% TODO: check inputs, transpose if necessary

n = size(vert,1);


if (ispc)
    D = inf(n);
    march = MESH.fastmarchmex('init', int32(tri-1), double(vert(:,1)), double(vert(:,2)), double(vert(:,3)));
    for i=1:n
        source = inf(n,1);
        source(i) = 0;
        D(:,i) = MESH.fastmarchmex('march', march, double(source));
    end
    MESH.fastmarchmex('deinit', march);
elseif (isunix)
    D = distmatrix(vert',tri');
else
    disp('No distance matrix calculated.');
    D=zeros(n);
end

% In case D is not symmetric
D = min(D,D');

end
