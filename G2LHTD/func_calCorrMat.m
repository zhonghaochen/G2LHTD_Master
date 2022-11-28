function [Rinv, R] = func_calCorrMat( imgVector )
[pixelCount, bandCount] = size( imgVector );
% ����ϵ������ R
R = 1/pixelCount * (imgVector' * imgVector);

% ��ȡ����ֵ�ԽǾ����������������
% �˴��൱��PCA�任��ȡǰPCs������ֵ����������
[eig_XL,eig_Z]=eig(R);
[Deig_Z,ind]=sort(diag(eig_Z),'descend');
D_eigXL=eig_XL(:,ind');

% �Զ�ȷ��ѡ������ɷָ���
rate = 0.9999;%�ò����ɵ�   
Sumva1 = rate * sum(Deig_Z); %���ܺ�0.99999������Сȡ������ֵ
T0=cumsum(Deig_Z);           % cumsumΪ�ۼӺ����������ۼ�  
ki=find(T0>Sumva1);   
PCs=ki(1);

Rinv=D_eigXL(:,1:PCs)*inv(diag(Deig_Z(1:PCs)))*D_eigXL(:,1:PCs)';
end

