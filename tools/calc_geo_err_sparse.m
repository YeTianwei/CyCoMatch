function errs = calc_geo_err_sparse(x,y, gt_matches, D)

% if (ndims(matches>1))
%    %TODO: 
% else

nm = length(y);
errs = zeros(nm,1);
for i=1:nm
    errs(i) = D( y(i), gt_matches(x(i)) );
end

errs = errs ./ max(max(D));

end