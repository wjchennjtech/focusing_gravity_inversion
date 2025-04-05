function [p] =Conjugradient2(A,dObs,Wd,Wm,ee,m0,Nmax,sigma,tolorence,m_low,m_max,mu)
global bounds;
global nx;
global ny;
global nz;
% mk=zeros(nx*ny*nz,1);
mk=m0;
invWm=diag(1./diag(Wm));
%添加SF函数
We=sparse(diag((mk.^2.+ee^2).^(-1/2).*((1+exp(-mk.^2./ee^2)).^(-1))));
invWe=sparse(diag((mk.^2.+ee^2).^(1/2).*(1+exp(-mk.^2./ee^2))));
%不添加SF函数
% We=sparse(diag((mk.^2.+ee^2).^(-1/2)));
% invWe=sparse(diag((mk.^2.+ee^2).^(1/2)));

dw=Wd*dObs;
% mw=Wm*We*(mk-m0);
mw=Wm*We*(mk);
Aw=Wd*A*invWm*invWe;

Q=invWe*((invWe))*(Aw')*Aw+mu*ee^2*eye(nx*ny*nz);
q=invWe*((invWe))*Aw'*dw;
f=Q*mw-q;
d=-f;
t=-(d'*f)/(d'*Q*d);
k=0;
while k<Nmax
    k=k+1;  
    mw=mw+t*d;
%     mk=invWm*invWe*mw+m0;
    mk=invWm*invWe*mw;
     if bounds==1
 mk = min(max(mk, m_low), m_max);
    end
  dd = norm(A * mk - dObs, 2);
    a(1,k)=dd;
    fg=figure(1);
    set(fg,'name','Residual','numbertitle','off')
    plot(1:k,a,'b-','linewidth',2)
    xlabel('Iterations');
    ylabel('Residuals');
    set(gca,'fontname','Times New Roman');
    if dd<tolorence
        break;
    end
    
    We=sparse(diag((mk.^2.+ee^2).^(-1/2).*((1+exp(-mk.^2./ee^2)).^(-1))));
    invWe=sparse(diag((mk.^2.+ee^2).^(1/2).*(1+exp(-mk.^2./ee^2))));
    %不添加SF函数
%     We=sparse(diag((mk.^2.+ee^2).^(-1/2)));
%     invWe=sparse(diag((mk.^2.+ee^2).^(1/2)));
    
%     mw=Wm*We*(mk-m0);
    mw=Wm*We*(mk);
    Aw=Wd*A*invWm*invWe;
    mu=sigma*mu;
    
    Q=invWe*((invWe))*(Aw')*Aw+mu*ee^2*eye(nx*ny*nz);
    q=invWe*((invWe))*Aw'*dw;
    fk=Q*mw-q;
    d=-fk+((fk'*fk)/(f'*f))*d;
    t=-(d'*fk)/(d'*Q*d);
    f=fk;    
end
res=a';
p=mk;
end



