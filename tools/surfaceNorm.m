function S = surfaceNorm_single(S)
X = [S.surface.X S.surface.Y S.surface.Z];

X = double(X);
n = size(X,1); 
normal.xm = mean(X);
X = X - repmat(normal.xm,n,1);
normal.xscale = sqrt(sum(sum(X.^2,2))/n);
X = X/normal.xscale;

S.surface.X = X(:,1);
S.surface.Y = X(:,2);
S.surface.Z = X(:,3);
S.surface.VERT = X;