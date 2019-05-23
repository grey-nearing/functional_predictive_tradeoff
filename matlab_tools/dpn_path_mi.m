function [Ixy,Hx,Hy] = dpn_path_te(X,Y,lag,Nb)

% get dimensions
X = X(:);
Y = Y(:);
[Nx,D] = size(Y); assert(D==1);
[Ny,D] = size(X); assert(D==1);
assert(Ny == Nx); N = Nx;

% deal with lags
Xw = window_average(X,lag); 
Yw = window_average(Y,lag);
assert(numel(Xw) == numel(Yw));
if numel(Xw)<1000; return; end;

% time shift
Xt = Xw(1:end-1,:); Xt = Xt(:);
Yt = Yw(1:end-1,:); Yt = Yt(:);
Ys = Yw(2:end,:);   Ys = Ys(:);
N = length(Ys);
assert(N==length(Xt));
assert(N==length(Yt));

% deal with grandmas
I = find(any(isnan([Xt,Yt,Ys]')));
Xt(I) = [];
Yt(I) = [];
Ys(I) = [];

% make sure we have dealt with missing values
assert(isempty(find(isnan(Xt),1,'first')));
assert(isempty(find(isnan(Yt),1,'first')));
assert(isempty(find(isnan(Ys),1,'first')));

% special case when targeting self
if max(abs(X-Y)==0)
 Yt = rand(size(Yt)) * (max(Yt)-min(Yt)) + min(Yt);
end

% bins
Bmin = min(Y)-1e-6; 
Bmax = max(Y)+1e-6; 
By = linspace(Bmin,Bmax,Nb);

Bmin = min(X)-1e-6;                
Bmax = max(X)+1e-6;                
Bx = linspace(Bmin,Bmax,Nb);

% do the actaul calculations
[Ixy,Hy,Hx] = mutual_info(Ys,Xt,By,Bx);

