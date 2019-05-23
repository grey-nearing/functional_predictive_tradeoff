function SI = Sobol(Y)

% dimensions
[Nsamples,Nparms] = size(Y);

% actual number of parameters (the frist run is not perturbed)
Nparms = Nparms-1;

% standardize (zero mean unit variance)
Y = (Y-repmat(mean(Y),Nsamples,1))./repmat(std(Y),Nsamples,1);

% non-perturbed results
Ysingle = Y(:,1);
fo = mean(Ysingle);
Do = var(Ysingle);
% Do = [mean(Y_single.^2)-fo^2]';

%% *** MC Variance Integrals **********************************************
D = zeros(Nparms,1);
for p = 1:Nparms
 for n = 1:Nsamples
  D(p) = D(p) + Y(n,1)*Y(n,p+1);
 end
end
D = D./Nsamples - fo^2;

% Total Effect
SI = 1-D./Do;

%DD_F = D(:,1)+D(:,3)./(Nsamples*2) - fo^2;
%DD_T = D(:,2)+D(:,4)./(Nsamples*2) - fo^2;
%FOSI = DD_F./Do;
%TSI  = 1 - DD_T./Do;

%% *** END FUNCTION *******************************************************
