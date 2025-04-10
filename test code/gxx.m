function [gxx] = gxx(xp,yp,zp,x1,x2,y1,y2,z1,z2,q)
%G is the gravitational constant
%Coordinates P1?(x1,y1,z1) and P2?(x2,y2,z2) represent the farthest diagonal corners of a rectangular prism
G=6.67e-11;
q1=q;
%Calculate the gravity anomaly of a finite-depth rectangular prism
gxx=0;
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
            gxx=(-1).^(s+d+f).*G.*q1.*(-atan(x(:,s).*y(:,d)./(x(:,s).*x(:,s)+r.*z(:,f)+(z(:,f).*z(:,f)))))+gxx;
        end
    end
end
gxx=-gxx;
end

