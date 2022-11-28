function [Rinv, R] = func_calCorrMat( imgVector )
[pixelCount, bandCount] = size( imgVector );
% 关联系数矩阵 R
R = 1/pixelCount * (imgVector' * imgVector);

% 求取特征值对角矩阵和特征向量矩阵
% 此处相当于PCA变换，取前PCs个特征值与特征向量
[eig_XL,eig_Z]=eig(R);
[Deig_Z,ind]=sort(diag(eig_Z),'descend');
D_eigXL=eig_XL(:,ind');

% 自动确定选择的主成分个数
rate = 0.9999;%该参数可调   
Sumva1 = rate * sum(Deig_Z); %按总和0.99999比例大小取舍特征值
T0=cumsum(Deig_Z);           % cumsum为累加函数，向下累加  
ki=find(T0>Sumva1);   
PCs=ki(1);

Rinv=D_eigXL(:,1:PCs)*inv(diag(Deig_Z(1:PCs)))*D_eigXL(:,1:PCs)';
end

