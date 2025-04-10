function [gx] = gx(xp,yp,zp,x1,x2,y1,y2,z1,z2,q)

G=6.67e-11;
q1=q;
gx=0;
x(:,1)=x1-xp;
x(:,2)=x2-xp;
y(:,1)=y1-yp;
y(:,2)=y2-yp;
z(:,1)=z1-zp;
z(:,2)=z2-zp;
for s=1:2
    for d=1:2
        for f=1:2
            %             a=xp-x(s);
            %             if a<1e-10
            %                 a=1e-10;
            %             else
            %                 a=a;
            %             end
            %             gx=(-1).^(s+d+f).*G.*q1.*((y(d)-yp).*log((z(f)-zp) + ((xp-x(s)).^2+(yp-y(d)).^2+(zp-z(f)).^2).^0.5) + (z(f)-zp).*log((y(d)-yp) + ((xp-x(s)).^2+(yp-y(d)).^2+(zp-z(f)).^2).^0.5)...
            %                            - (x(s)-xp).*atan((z(f)-zp).*(y(d)-yp)./a./((xp-x(s)).^2+(yp-y(d)).^2+(zp-z(f)).^2).^0.5))+gx ;
            r=sqrt(x(:,s).*x(:,s)+y(:,d).*y(:,d)+z(:,f).*z(:,f));
            gx=(-1).^(s+d+f).*G.*q1.*((y(:,d).*log(z(:,f) + r) + z(:,f).*log(y(:,d) + r) ...
                           - x(:,s).*atan(z(:,f).*y(:,d)./x(:,s)./r)))+gx;
        end
    end
end
gx = fillmissing(gx,'previous');
end

