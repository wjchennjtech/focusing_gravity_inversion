function [p] =Conjugradient3(A,dObs,Wd,Wm,ee,m0,Nmax,sigma,tolorence,m_low,m_max,mu)
global bounds nx ny nz;
% mk=zeros(nx*ny*nz,1);
mk=m0;
invWm=diag(1./diag(Wm));
%add SF function
We=sparse(diag((mk.^2.+ee^2).^(-1/2).*((1+exp(-mk.^2./ee^2)).^(-1))));
invWe=sparse(diag((mk.^2.+ee^2).^(1/2).*(1+exp(-mk.^2./ee^2))));

% We=sparse(diag((mk.^2.+ee^2).^(-1/2)));
% invWe=sparse(diag((mk.^2.+ee^2).^(1/2)));

% Calculate Aw, Q, q
dw = Wd * dObs;
mw=Wm*We*(mk);
Aw = Wd * A * invWm * invWe;
Q=invWe*((invWe))*(Aw')*Aw+mu*ee^2*eye(nx*ny*nz);
q=invWe*((invWe))*Aw'*dw;

% Compute the initial gradient and search direction
f = Q * (Wm * We * mk) - q;
M = diag(Q);  % Preconditioning Matrix
Minv = sparse(diag(1 ./ (M + 1e-6))); 
d = -Minv * f;  % Preconditioned Conjugate Directions
% d=-f;
t = -(d' * f) / (d' * Q * d);

k = 0;
resHist = zeros(Nmax, 1);  % record residuals

while k<Nmax
    k=k+1;
    mw=mw+t*d;
    %     mk=invWm*invWe*mw+m0;
    mk=invWm*invWe*mw;
    if bounds == 1
        mk = min(max(mk, m_low), m_max);
    end
    % calculate residuals
    res = norm(A * mk - dObs, 2);  %L2 norm calculation
    resHist(k) = res;
    % visualize the residual descent process
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
   
    %     We=sparse(diag((mk.^2.+ee^2).^(-1/2)));
    %     invWe=sparse(diag((mk.^2.+ee^2).^(1/2)));
    
    %     mw=Wm*We*(mk-m0);
    mw=Wm*We*(mk);
    mu=sigma*mu;
    Aw = Wd * A * invWm * invWe;
   Q=invWe*((invWe))*(Aw')*Aw+mu*ee^2*eye(nx*ny*nz);
    q=invWe*((invWe))*Aw'*dw;
    fk=Q*mw-q;
    % ¼ÆËã ¦Â_k
    beta = (fk' * fk) / (f' * f);
    % update the search direction
    d = -Minv * fk + beta * d;
        d=-fk+((fk'*fk)/(f'*f))*d;
    t=-(d'*fk)/(d'*Q*d);
%     **adaptive step size**
    if k > 1
        if resHist(k) > resHist(k-1) 
            t = t * 0.5;
        elseif resHist(k) < resHist(k-1) * 0.9 
            t = t * 1.1;
        end
    end
    f=fk;
end
p=mk;
end
