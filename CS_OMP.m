function [x,Pos_x] = CS_OMP(y,A,t)
%CS_OMP 此处显示有关此函数的摘要
%   y = A*x，其中A: R*K，y:R*1，x:K*1
%   R个接收机，t个未知用户，t为稀疏度
%   现在已知y和A，求x
    [y_rows,y_columns] = size(y);  
    if y_rows < y_columns
        y = y';   % y should be a column vector
    end
    [R,K] = size(A);  % A: sensing matrix
    x = zeros(K,1);  % x is a column vector to be predicted
    At = zeros(R,t);  % At用来存储迭代过程中A被选择的列
    Pos_x = zeros(1,t);  % Pos_x用来存储A中被选择的列的序号
    r_n = y;  % initial residual is y
    for ii = 1:t
        product = A'*r_n;  % 计算传感矩阵A各列与残差的内积，K*1
        [val,pos] = max(abs(product));   % 寻找最大内积绝对值，即与残差最相关的列
        At(:,ii) = A(:,pos);
        Pos_x(ii) = pos;
        A(:,pos) = zeros(R,1);  %清零A的这一列，可不要
        x_ls = (At(:,1:ii)'*At(:,1:ii))^(-1)*At(:,1:ii)'*y;  % 求最小二乘解
        r_n = y - At(:,1:ii)*x_ls;  % 更新残差
    end
    x(Pos_x) = x_ls;  % 恢复出来的x
end

