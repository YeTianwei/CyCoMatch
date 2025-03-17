
clear;clc;close all;

addpath(genpath('./'));

% load shapes

Spath = './data/';
finf = dir(fullfile(Spath, '*.ply'));

numShapes = length(finf);
shapes = cell(numShapes, 1);
D = cell(numShapes, 1);

k = 500;
opts.shot_num_bins = 10;
opts.shot_radius = 5;

for i = 1:numShapes
    file_path = fullfile(Spath, finf(i).name);
    shapes{i} = MESH.preprocess(file_path, 'IfComputeLB', true, 'numEigs', k,'IfFindNeigh',true,'IfFindEdge',true,'IfComputeGeoDist',true,'IfComputeNormals',true);
    shapes{i} = surfaceNorm(shapes{i});
    shapes{i}.area = sum(calc_tri_areas(shapes{i}.surface));
    shapes{i}.SHOT = calc_shot(shapes{i}.surface.VERT', shapes{i}.surface.TRIV', 1:shapes{i}.nv, opts.shot_num_bins, opts.shot_radius*sqrt(shapes{i}.area)/100, 3)';
    D{i} = shapes{i}.Gamma;    
end

% compute T0 as initialization

T0 = cell(numShapes);
for i = 1:numShapes
    for j = 1:numShapes
        T0{i,j} = knnsearch(gpuArray(shapes{j}.SHOT), gpuArray(shapes{i}.SHOT)); 
        T0{i,j} = gather(T0{i,j});
        if i == j,T0{i,j} = []; end
    end
end

G = 1 - cellfun(@isempty, T0);

for i = 1:numShapes, T0{i,i} = (1:shapes{i}.nv)'; end

% main

[C, T] = Multi_CCB(shapes, T0, D, G);

% plot results
figure;
for idx = 1:numShapes
    subplot(2,3,idx);
    rgb = coord2rgb(shapes{idx}.surface.VERT);
    mplot_mesh_rgb(shapes{idx}.surface.VERT,shapes{idx}.surface.TRIV,rgb(T{1,idx}));
end



