function [p] =Conjugradient3(A,dObs,Wd,Wm,ee,m0,Nmax,sigma,tolorence,m_low,m_max,mu)
global bounds nx ny nz;
% mk=zeros(nx*ny*nz,1);
mk=m0;
invWm=diag(1./diag(Wm));
%���SF����
We=sparse(diag((mk.^2.+ee^2).^(-1/2).*((1+exp(-mk.^2./ee^2)).^(-1))));
invWe=sparse(diag((mk.^2.+ee^2).^(1/2).*(1+exp(-mk.^2./ee^2))));
%�����SF����
% We=sparse(diag((mk.^2.+ee^2).^(-1/2)));
% invWe=sparse(diag((mk.^2.+ee^2).^(1/2)));
% ���� Aw, Q, q
dw = Wd * dObs;
mw=Wm*We*(mk);
Aw = Wd * A * invWm * invWe;
Q=invWe*((invWe))*(Aw')*Aw+mu*ee^2*eye(nx*ny*nz);
q=invWe*((invWe))*Aw'*dw;

% �����ʼ�ݶȺ���������
f = Q * (Wm * We * mk) - q;
M = diag(Q);  % Ԥ��������
Minv = sparse(diag(1 ./ (M + 1e-6)));  % Ԥ������
d = -Minv * f;  % Ԥ���������
% d=-f;
t = -(d' * f) / (d' * Q * d);

k = 0;
resHist = zeros(Nmax, 1);  % ��¼�в�

while k<Nmax
    k=k+1;
    mw=mw+t*d;
    %     mk=invWm*invWe*mw+m0;
    mk=invWm*invWe*mw;
    if bounds == 1
        mk = min(max(mk, m_low), m_max);
    end
    % ����в�
    res = norm(A * mk - dObs, 2);  % L2 ��������
    resHist(k) = res;
    % ���ӻ��в��½�����
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
    %�����SF����
    %     We=sparse(diag((mk.^2.+ee^2).^(-1/2)));
    %     invWe=sparse(diag((mk.^2.+ee^2).^(1/2)));
    
    %     mw=Wm*We*(mk-m0);
    mw=Wm*We*(mk);
    mu=sigma*mu;
    Aw = Wd * A * invWm * invWe;
   Q=invWe*((invWe))*(Aw')*Aw+mu*ee^2*eye(nx*ny*nz);
    q=invWe*((invWe))*Aw'*dw;
    fk=Q*mw-q;
    % ���� ��_k
    beta = (fk' * fk) / (f' * f);
    % ������������
    d = -Minv * fk + beta * d;
        d=-fk+((fk'*fk)/(f'*f))*d;
    t=-(d'*fk)/(d'*Q*d);
%     **����Ӧ��������**
    if k > 1
        if resHist(k) > resHist(k-1)  % ����в����ӣ���������
            t = t * 0.5;
        elseif resHist(k) < resHist(k-1) * 0.9  % ����в��½��죬��������
            t = t * 1.1;
        end
    end
    f=fk;
end
p=mk;

% %%
% % ��ʼ������
% mk = m0;  
% invWm = diag(1 ./ diag(Wm));  % ���� Wm ����ԽǾ���
% 
% % ����Ȩ�� We ��������� invWe
% We = sparse(diag((mk.^2 + ee^2).^(-1/2) .* ((1 + exp(-mk.^2 ./ ee^2)).^(-1))));
% invWe = sparse(diag((mk.^2 + ee^2).^(1/2) .* (1 + exp(-mk.^2 ./ ee^2))));
% 
% % �������ݲв�
% dw = Wd * dObs;
% mw = Wm * We * (mk);  % �任����Ȩ������
% Aw = Wd * A * invWm * invWe;  % ��Ȩ�����ռ��е�ϵͳ����
% 
% % ��ʼ���������� d ���ݶ� f
% Q = invWe * (invWe' * (Aw' * Aw)) + mu * ee^2 * eye(nx * ny * nz);
% q = invWe * (invWe' * Aw' * dw);
% f = Q * mw - q;
% d = -f;
% t = -(d' * f) / (d' * Q * d);  % ���㲽��
% 
% % ��������
% k = 0;
% resHist = zeros(maxIter, 1);  % ��¼�в���ʷ
% 
% while k < maxIter
%     k = k + 1;
% 
%     % ���¼�Ȩ������Ĳ��� mw
%     mw = mw + t * d;
%     
%     % ���任�����Բ�����
%     mk = invWm * invWe * mw;  
% 
%     % Լ�����±߽�
%     if bounds == 1
%         mk = min(max(mk, m_low), m_max);
%     end
%     
%     % ����в�
%     res = norm(A * mk - dObs, 2);  % L2 ��������
%     resHist(k) = res;
%     
%     % ���ӻ��в��½�����
%     figure(1);
%     set(gcf, 'name', 'Residual', 'numbertitle', 'off');
%     plot(1:k, resHist(1:k), 'b-', 'linewidth', 2);
%     xlabel('Iterations');
%     ylabel('Residuals');
%     set(gca, 'FontName', 'Times New Roman');
%     
%     % ��ֹ����
%     if res < tol
%         break;
%     end
% 
%     % ���� We
%     We = sparse(diag((mk.^2 + ee^2).^(-1/2) .* ((1 + exp(-mk.^2 ./ ee^2)).^(-1))));
%     invWe = sparse(diag((mk.^2 + ee^2).^(1/2) .* (1 + exp(-mk.^2 ./ ee^2))));
%     
%     % ���� Q �� q
%     mw = Wm * We * (mk);
%     Aw = Wd * A * invWm * invWe;
%     
%     mu = sigma * mu;  % ����Ӧ�������򻯲���
%     
%     Q = invWe * (invWe' * (Aw' * Aw)) + mu * ee^2 * eye(nx * ny * nz);
%     q = invWe * (invWe' * Aw' * dw);
%     
%     % �����µ��ݶ� fk
%     fk = Q * mw - q;
%     
%     % ������������
%     beta = (fk' * fk) / (f' * f);
%     d = -fk + beta * d;
%     
%     % ���㲽�� t
%     t = -(d' * fk) / (d' * Q * d);
% 
%     % **����Ӧ��������**
%     if k > 1
%         if resHist(k) > resHist(k-1)  % ����в����ӣ���������
%             t = t * 0.5;
%         elseif resHist(k) < resHist(k-1) * 0.9  % ����в��½��죬��������
%             t = t * 1.1;
%         end
%     end
%     
%     f = fk;  % �����ݶ�
% end
% 
% % ������շ��ݽ��
% p = mk;
end
