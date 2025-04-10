function [gzz] = gzz(xp,yp,zp,x1,x2,y1,y2,z1,z2,q)

G=6.67e-11;
q1=q;
gzz=0;
x(:,1)=x1-xp;
x(:,2)=x2-xp;
y(:,1)=y1-yp;
y(:,2)=y2-yp;
z(:,1)=z1-zp;
z(:,2)=z2-zp;
for s=1:2
    for d=1:2
        for f=1:2
            r=sqrt(x(:,s).*x(:,s)+y(:,d).*y(:,d)+z(:,f).*z(:,f));
            gzz=(-1).^(s+d+f).*G.*q1.*(-atan(x(:,s).*y(:,d)./-z(:,f)./r))+gzz;
        end
    end
end
gzz=-gzz;
end
