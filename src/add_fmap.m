function C = add_fmap(C1,C2,alpha)

n = length(C1);
C = cell(n);

dim = 500;

for i = 1:n
    for j = 1:n
        
        if i == j,C{i,j} = C1{i,j}; continue;  end
        
        [rowsC2, colsC2] = size(C2{i,j});
        
        padRows = dim - rowsC2;
        padCols = dim - colsC2;

        if padRows < 0 || padCols < 0
            error('C2{%d,%d} exceeds the specified dimension %d', i, j, dim);
        end
        
        C2_padded = padarray(C2{i,j}, [0, padCols], 0, 'post'); 
        C2_padded = padarray(C2_padded, [padRows, 0], 0, 'post');  
        
        [rowsC1, colsC1] = size(C1{i,j});
        if rowsC1 ~= dim || colsC1 ~= dim
            error('C1{%d,%d} must be of size %dx%d', i, j, dim, dim);
        end
        
        C{i,j} = alpha * C1{i,j} + (1 - alpha) * C2_padded;
        
    end
end
end