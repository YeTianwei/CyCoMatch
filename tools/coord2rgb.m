function rgb=coord2rgb(V)
% 将输入矩阵 V 中的值进行归一化到 [0, 1] 的范围

min_value=min(V,[],1);
max_value=max(V,[],1);

rgb=(V-min_value)./(max_value-min_value);

