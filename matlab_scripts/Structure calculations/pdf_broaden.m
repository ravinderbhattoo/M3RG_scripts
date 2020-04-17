function [Q,SQ]=pdf_broaden(xbins,pdf,qrange,dq,Rmax)

Q=((floor(qrange(1)/dq):floor(qrange(2)/dq))-0.5)*dq;
Q=qrange(1):dq:qrange(2)-dq;
    fwhm=5.437/Rmax;
    [~,Iend]=min(abs(xbins-Rmax));

SQ=zeros(1,length(Q));
for i=1:length(pdf)
    for j=1:length(Q) 
        SQ(j)=SQ(j)+exp(-((i-j)*dq)^2/(2*(fwhm/2.35482)^2))*pdf(i)*(1/sqrt(2*pi*(fwhm/2.35482)));       
    end
end
