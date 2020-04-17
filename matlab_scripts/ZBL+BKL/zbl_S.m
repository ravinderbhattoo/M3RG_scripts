function E=zbl_S(r,Z1,Z2,rin,rout)

dd=0.001;

rfit=linspace(rout-dd,rout+dd,4);
pp=polyfit(rfit,zbl(rfit,Z1,Z2),2);

% figure;
% plot(rfit,zbl(rfit,Z1,Z2),'.',rfit,polyval(pp,rfit))
% figure;

pp1=polyder(pp);
pp2=polyder(pp1);

Er=zbl(rout,Z1,Z2);
dEr=polyval(pp1,rout);
ddEr=polyval(pp2,rout);

A=(-3*dEr+(rout-rin)*ddEr)/(rout-rin)^2;
B=(2*dEr-(rout-rin)*ddEr)/(rout-rin)^3;
C=-Er+0.5*(rout-rin)*dEr-(rout-rin)^2*ddEr/12;

S=zeros(size(r));
S(r<=rin)=C;
S(r>rin)=A/3*(r(r>rin)-rin).^3+B/4*(r(r>rin)-rin).^4+C;
E=zbl(r,Z1,Z2)+S;