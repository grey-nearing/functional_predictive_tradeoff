function Xw = window_average(X,W)

% get dimensions
[N,D] = size(X);
assert(N > W*100);
assert(D == 1);

% bin centers
Ic = floor(W/2)+1:W:N-floor(W/2);
Nc = length(Ic);

% init storage
Xw = zeros(Nc,1)./0;

% loop through dimensions of input data
for c = 1:Nc
 sdex = Ic(c) - floor(W/2);
 edex = Ic(c) + floor(W/2);
% [sdex,edex,length(sdex:edex)]
 if length(find(isnan(X(sdex:edex)))) < 0.1*W
  Xw(c) = nanmean(X(sdex:edex));
 end
end

