clear
clc

yeo7 = xlsread('Schaefer2018_200Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv');
yeo7(:,1:2)=[];

yeo17 = xlsread('Schaefer2018_200Parcels_17Networks_order_FSLMNI152_1mm.Centroid_RAS.csv');
yeo17(:,1:2)=[];

for i = 1:200
    idx(:,1) = (yeo7(:,1) == yeo17(i,1));
    idx(:,2) = (yeo7(:,2) == yeo17(i,2));
    idx(:,3) = (yeo7(:,3) == yeo17(i,3));
    Schaefer200_7to17(i,1) = find(sum(idx,2)==3);

    idx(:,1) = (yeo17(:,1) == yeo7(i,1));
    idx(:,2) = (yeo17(:,2) == yeo7(i,2));
    idx(:,3) = (yeo17(:,3) == yeo7(i,3));
    Schaefer200_17to7(i,1) = find(sum(idx,2)==3);
end

length(unique(Schaefer200_7to17))
length(unique(Schaefer200_17to7))

cd('F:\Cui_Lab\Projects\Connectional_Variability_Axis\functions')
save('Schaefer200_7to17.mat','Schaefer200_7to17')
save('Schaefer200_17to7.mat','Schaefer200_17to7')