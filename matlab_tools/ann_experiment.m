function [model,results,Yhat] = ann_experiment(X,Y,I,J,trnfctn,Nepochs,verbose)

% make sure no missing values
assert(isempty(find(isnan(X(:)))));
assert(isempty(find(isnan(Y(:)))));
assert(isempty(find(X(:)<-9990)));
assert(isempty(find(Y(:)<-9990)));

% train
if verbose
 fprintf('\nTraining: Ntrn = %d - Nepochs %d ... ',length(I),Nepochs); tic;
end
model = train_ann(X,Y,I,trnfctn,Nepochs);
if verbose
 t = toc; fprintf('finished - time = %f \n\n',t);
end

% test
if verbose
 fprintf('Predicting: Ndata = %d ... ',length(Y)); tic;
end
[results,Yhat] = pred_ann(X,Y,model,I,J);
if verbose
 t = toc; fprintf('finished - time = %f \n\n',t);
end

% screen report
if verbose
 fprintf('Ntrn = %d    Dx = %d \n',length(I),size(X,2));
 fprintf(' - Info (1%%) : %f, %f, %f \n',squeeze(results.info(:,1)));
 fprintf(' - Info (2%%) : %f, %f, %f \n',squeeze(results.info(:,2)));
 fprintf(' - Info (5%%) : %f, %f, %f \n',squeeze(results.info(:,3)));
 fprintf(' - Info (10%%): %f, %f, %f \n',squeeze(results.info(:,4)));
 fprintf(' - RMSE      : %f, %f, %f \n' ,squeeze(results.rmse));
 fprintf(' - CORR      : %f, %f, %f \n' ,squeeze(results.corr));
 fprintf('\n');
end

