clear all;
load Sandiego;

addpath(genpath('./EMAP'));

%% 数据处理(PCA)
[row,col,bands] = size(data);
data_scaled = hyperNormalize(data);

data_scaled_2d = hyperConvert2d(data_scaled);
num_of_pca = 3;

data_pca = pca(data_scaled_2d',num_of_pca);
data_pca_3d = hyperConvert3d(data_pca',row,col,num_of_pca);


%% EMAP
tic
data_EMAP = EMAP(data_pca_3d,'data_EMAP',false, false,'a', [100,500,1000,5000],'d',[10,25,50,100],'i', [0.2,0.3,0.4,0.5],'s',[20,30,40,50]);
toc

data_EMAP = double(data_EMAP);
data_EMAP = hyperNormalize(data_EMAP);
%% 获取目标光谱
target = target_emap;
%% D2CEM 
win_out=25;
win_in=7;
tic
results = D2CEM(data_EMAP,target,win_out,win_in);
toc
results = hyperNormalize(results);

%% 结果展示
[FPR_1_our,TPR_1_our,thre_1_our,auc_1_our] = myPlot3DROC(map,results);
figure;imagesc(map);colormap('gray'); axis image;  title('groundtruth') ; set(gca,'xtick',[]); set(gca,'ytick',[]);set(gca,'Position',[0 0 1 1]);set (gcf,'Position',[0,0,512,512]);
figure;imagesc(results); colormap('gray');axis image;  title('results') ; set(gca,'xtick',[]); set(gca,'ytick',[]);set(gca,'Position',[0 0 1 1]);set (gcf,'Position',[0,0,512,512]);


