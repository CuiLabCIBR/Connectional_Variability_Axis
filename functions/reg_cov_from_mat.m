function data_reg_cov = reg_cov_from_mat(data_mat,cov_mat,data_mask)
% This function is used to regress the covariates from the connectivity matrix.
% input:
%       data_mat: N*N connectivity matrix, N is the number of nodes.
%       cov_mat: N*N*M covariates matrix, M is the number of covariates.
%       data_mask: N*N matrix, analyze the data within the mask.
% output:
%       data_reg_cov: N*N connectivity matrix, regressed the covariates matrix.


N = length(data_mat);
[~,~,M] = size(cov_mat);

if ~exist('data_mask','var')
    data_mask = ones(N,N);
end

for i = 1:M
    cov_vec(:,i) = mat2vec(cov_mat(:,:,i));
end

data_vec = mat2vec(data_mat); % extract the lower triangle elements as a vector
mask_vec = mat2vec(data_mask); % extract the lower triangle elements as a vector

x = cov_vec(mask_vec>0,:);
x = bsxfun(@minus, x, mean(x));

% Regression
y = data_vec(mask_vec>0); 
[b,dev,stats] = glmfit(x,y); %% regression
con = b(1); %constsant term
resid = stats.resid; % residuals
z = resid+con; %% add constant term

data_vec(mask_vec>0) = z;
data_reg_cov = squareform(data_vec); % convert vector to matrix
