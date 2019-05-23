function [Iyx,Hy] = transfer_info(Y,X,Bx,By)

[N,D] = size(X); assert(N == size(Y,1)); assert(D==1);

assert(isempty(find(isnan(Y))));
assert(isempty(find(isnan(X))));

assert(isempty(find(Y<=-9990)));
assert(isempty(find(X<=-9990)));

%Bmin = min(X)-1e-6; 
%Bmax = max(X)+1e-6; 
%Bx = Bmin:dx:Bmax;
%
%Bmin = min(Y)-1e-6; 
%Bmax = max(Y)+1e-6; 
%By = Bmin:dy:Bmax;

[Pyx,Py,Px] = hist2(Y,X,By,Bx);

Pyx  = Psx(:);
Py   = Ps(:);
Px   = Px(:);

if abs(sum(Py)-1)>1/N^2; error(' ');end;%('Ps does not sum to 1'); end;
if abs(sum(Px)-1)>1/N^2; error(' ');end;%('Px does not sum to 1'); end;

if abs(sum(Pyx)-1)>1/N^2; error(' ');end;%('Psx does not sum to 1'); end;

Hy = -Py(Py>0)'*log(Py(Py>0));
Hx = -Px(Px>0)'*log(Px(Px>0));
Hyx = -Pyx(Pyx>0)'*log(Pyx(Pyx>0));

Iyx = Hy+Hx-Hyx;






