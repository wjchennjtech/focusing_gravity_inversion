function [p] =Conjugradient3(A,dObs,Wd,Wm,ee,m0,Nmax,sigma,tolorence,m_low,m_max,mu)
global bounds nx ny nz;
% mk=zeros(nx*ny*nz,1);
mk=m0;
invWm=diag(1./diag(Wm));
%添加SF函数
We=sparse(diag((mk.^2.+ee^2).^(-1/2).*((1+exp(-mk.^2./ee^2)).^(-1))));
invWe=sparse(diag((mk.^2.+ee^2).^(1/2).*(1+exp(-mk.^2./ee^2))));
%不添加SF函数
% We=sparse(diag((mk.^2.+ee^2).^(-1/2)));
% invWe=sparse(diag((mk.^2.+ee^2).^(1/2)));
% 计算 Aw, Q, q
dw = Wd * dObs;
mw=Wm*We*(mk);
Aw = Wd * A * invWm * invWe;
Q=invWe*((invWe))*(Aw')*Aw+mu*ee^2*eye(nx*ny*nz);
q=invWe*((invWe))*Aw'*dw;

% 计算初始梯度和搜索方向
f = Q * (Wm * We * mk) - q;
M = diag(Q);  % 预条件矩阵
Minv = sparse(diag(1 ./ (M + 1e-6)));  % 预条件逆
d = -Minv * f;  % 预条件共轭方向
% d=-f;
t = -(d' * f) / (d' * Q * d);

k = 0;
resHist = zeros(Nmax, 1);  % 记录残差

while k<Nmax
    k=k+1;
    mw=mw+t*d;
    %     mk=invWm*invWe*mw+m0;
    mk=invWm*invWe*mw;
    if bounds == 1
        mk = min(max(mk, m_low), m_max);
    end
    % 计算残差
    res = norm(A * mk - dObs, 2);  % L2 范数计算
    resHist(k) = res;
    % 可视化残差下降过程
    figure(1);
    set(gcf, 'name', 'Residual', 'numbertitle', 'off');
    plot(1:k, resHist(1:k), 'b-', 'linewidth', 2);
    xlabel('Iterations');
    ylabel('Residuals');
    set(gca, 'FontName', 'Times New Roman');
    if res<tolorence
        break;
    end
    We=sparse(diag((mk.^2.+ee^2).^(-1/2).*((1+exp(-mk.^2./ee^2)).^(-1))));
    invWe=sparse(diag((mk.^2.+ee^2).^(1/2).*(1+exp(-mk.^2./ee^2))));
    %不添加SF函数
    %     We=sparse(diag((mk.^2.+ee^2).^(-1/2)));
    %     invWe=sparse(diag((mk.^2.+ee^2).^(1/2)));
    
    %     mw=Wm*We*(mk-m0);
    mw=Wm*We*(mk);
    mu=sigma*mu;
    Aw = Wd * A * invWm * invWe;
   Q=invWe*((invWe))*(Aw')*Aw+mu*ee^2*eye(nx*ny*nz);
    q=invWe*((invWe))*Aw'*dw;
    fk=Q*mw-q;
    % 计算 β_k
    beta = (fk' * fk) / (f' * f);
    % 更新搜索方向
    d = -Minv * fk + beta * d;
        d=-fk+((fk'*fk)/(f'*f))*d;
    t=-(d'*fk)/(d'*Q*d);
%     **自适应步长调整**
    if k > 1
        if resHist(k) > resHist(k-1)  % 如果残差增加，步长减半
            t = t * 0.5;
        elseif resHist(k) < resHist(k-1) * 0.9  % 如果残差下降快，步长增大
            t = t * 1.1;
        end
    end
    f=fk;
end
p=mk;

% %%
% % 初始化参数
% mk = m0;  
% invWm = diag(1 ./ diag(Wm));  % 计算 Wm 的逆对角矩阵
% 
% % 计算权重 We 和其逆矩阵 invWe
% We = sparse(diag((mk.^2 + ee^2).^(-1/2) .* ((1 + exp(-mk.^2 ./ ee^2)).^(-1))));
% invWe = sparse(diag((mk.^2 + ee^2).^(1/2) .* (1 + exp(-mk.^2 ./ ee^2))));
% 
% % 计算数据残差
% dw = Wd * dObs;
% mw = Wm * We * (mk);  % 变换到加权参数域
% Aw = Wd * A * invWm * invWe;  % 加权参数空间中的系统矩阵
% 
% % 初始化搜索方向 d 和梯度 f
% Q = invWe * (invWe' * (Aw' * Aw)) + mu * ee^2 * eye(nx * ny * nz);
% q = invWe * (invWe' * Aw' * dw);
% f = Q * mw - q;
% d = -f;
% t = -(d' * f) / (d' * Q * d);  % 计算步长
% 
% % 迭代参数
% k = 0;
% resHist = zeros(maxIter, 1);  % 记录残差历史
% 
% while k < maxIter
%     k = k + 1;
% 
%     % 更新加权参数域的参数 mw
%     mw = mw + t * d;
%     
%     % 反变换回物性参数域
%     mk = invWm * invWe * mw;  
% 
%     % 约束上下边界
%     if bounds == 1
%         mk = min(max(mk, m_low), m_max);
%     end
%     
%     % 计算残差
%     res = norm(A * mk - dObs, 2);  % L2 范数计算
%     resHist(k) = res;
%     
%     % 可视化残差下降过程
%     figure(1);
%     set(gcf, 'name', 'Residual', 'numbertitle', 'off');
%     plot(1:k, resHist(1:k), 'b-', 'linewidth', 2);
%     xlabel('Iterations');
%     ylabel('Residuals');
%     set(gca, 'FontName', 'Times New Roman');
%     
%     % 终止条件
%     if res < tol
%         break;
%     end
% 
%     % 更新 We
%     We = sparse(diag((mk.^2 + ee^2).^(-1/2) .* ((1 + exp(-mk.^2 ./ ee^2)).^(-1))));
%     invWe = sparse(diag((mk.^2 + ee^2).^(1/2) .* (1 + exp(-mk.^2 ./ ee^2))));
%     
%     % 计算 Q 和 q
%     mw = Wm * We * (mk);
%     Aw = Wd * A * invWm * invWe;
%     
%     mu = sigma * mu;  % 自适应调整正则化参数
%     
%     Q = invWe * (invWe' * (Aw' * Aw)) + mu * ee^2 * eye(nx * ny * nz);
%     q = invWe * (invWe' * Aw' * dw);
%     
%     % 计算新的梯度 fk
%     fk = Q * mw - q;
%     
%     % 计算搜索方向
%     beta = (fk' * fk) / (f' * f);
%     d = -fk + beta * d;
%     
%     % 计算步长 t
%     t = -(d' * fk) / (d' * Q * d);
% 
%     % **自适应步长调整**
%     if k > 1
%         if resHist(k) > resHist(k-1)  % 如果残差增加，步长减半
%             t = t * 0.5;
%         elseif resHist(k) < resHist(k-1) * 0.9  % 如果残差下降快，步长增大
%             t = t * 1.1;
%         end
%     end
%     
%     f = fk;  % 更新梯度
% end
% 
% % 输出最终反演结果
% p = mk;
end
