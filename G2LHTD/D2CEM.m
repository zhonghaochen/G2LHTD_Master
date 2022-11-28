function ResImg=D2CEM(hsi,targetSpec,win_out,win_in)
%% 局部四方向
[rows,cols,bands] = size(hsi);
%%
% 边界扩展
t = fix(win_out/2);
t1 = fix(win_in/2);

DataTest = zeros(3*rows, 3*cols, bands);
DataTest(rows+1:2*rows, cols+1:2*cols, :) = hsi;
DataTest(rows+1:2*rows, 1:cols, :) = hsi(:,cols:-1:1, :);
DataTest(rows+1:2*rows, 2*cols+1:3*cols, :) = hsi(:, cols:-1:1, :);
DataTest(1:rows, :, :) = DataTest(2*rows:-1:(rows+1), :, :);
DataTest(2*rows+1:3*rows, :, :) = DataTest(2*rows:-1:(rows+1), :, :);
target = targetSpec;
RDet=zeros(rows,cols);
for i = 1+rows: 2*rows      % 行
    for j = 1+cols: 2*cols  % 列
        aa = [i,i-1,i+1,i,i,i-1,i+1,i+1,i-1];
        bb = [j,j,j,j-1,j+1,j-1,j+1,j-1,j+1];
        ResCEM = zeros(1,9);
        for cc = 1:9
            block = DataTest(aa(cc)-t:aa(cc)+t,bb(cc)-t:bb(cc)+t,:);
            block(t-t1+1:t+t1+1, t-t1+1:t+t1+1, :) = NaN;
            block = reshape(block, win_out*win_out, bands);
            block(isnan(block(:, 1)), :) = [];
            imgVector = block;  % num_sam x num_dim
            CenPix= squeeze(DataTest(aa(cc), bb(cc), :));   % dim x 1 
            RInv = func_calCorrMat(imgVector); % correlation matrix
            w = RInv*targetSpec/(targetSpec'*RInv*targetSpec);
            ResCEM(cc)=w'*CenPix;
        end
        % 皮尔森相关系数
        p1 = fun_pear(squeeze(DataTest(i,j,:)),target);    

        px2 = fun_pear(squeeze(DataTest(i-1,j,:)),squeeze(DataTest(i,j,:)));
        px3 = fun_pear(squeeze(DataTest(i+1,j,:)),squeeze(DataTest(i,j,:)));
        px4 = fun_pear(squeeze(DataTest(i,j-1,:)),squeeze(DataTest(i,j,:)));
        px5 = fun_pear(squeeze(DataTest(i,j+1,:)),squeeze(DataTest(i,j,:)));
        px6 = fun_pear(squeeze(DataTest(i-1,j-1,:)),squeeze(DataTest(i,j,:)));
        px7 = fun_pear(squeeze(DataTest(i+1,j+1,:)),squeeze(DataTest(i,j,:)));
        px8 = fun_pear(squeeze(DataTest(i+1,j-1,:)),squeeze(DataTest(i,j,:)));
        px9 = fun_pear(squeeze(DataTest(i-1,j+1,:)),squeeze(DataTest(i,j,:)));
              
        % SAM值
        s1 = hyperSam(squeeze(DataTest(i,j,:)),target);

        % 8邻域ANFA       
        neighbour = px2*ResCEM(2)+px3*ResCEM(3)+px4*ResCEM(4)+px5*ResCEM(5)+px6*ResCEM(6)+px7*ResCEM(7)+px8*ResCEM(8)+px9*ResCEM(9);
        results = neighbour+ p1/s1*ResCEM(1);
        RDet(i-rows,j-cols) = results;
    end
end
ResImg = RDet;
end
function coeff = fun_pear(X,Y)  
    fenzi = sum(X.* Y) - (sum(X) * sum(Y)) / length(X);  
    fenmu = sqrt((sum(X.^2) - sum(X)^2 / length(X)) * (sum(Y.^2) - sum(Y)^2 / length(X)));  
    coeff = fenzi / fenmu; 
end

