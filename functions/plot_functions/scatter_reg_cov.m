function scatter_reg_cov(y1,y2,cov)
% 
num = length(y1);
[r,p] = partialcorr(y1,y2,cov);

[b1,~,residual_y1] = regress(y1,[ones(num,1),cov]);
[b2,~,residual_y2] = regress(y2,[ones(num,1),cov]);

% figure
% scatter(b1(1)+residual_y1,b2(1)+residual_y2)
% corr(residual_y1,residual_y2)

X = b1(1)+residual_y1;
Y = b2(1)+residual_y2;
mdl = fitlm(X,Y);
plot(mdl)