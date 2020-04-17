function [Q,SQ]=SQfz(xbins,pdf,pho,qrange,dq,use_window,varargin)

% 10/05/2015, Bu

Q=((floor(qrange(1)/dq):floor(qrange(2)/dq))-0.5)*dq;

if use_window==1
    Rmax=varargin{1};
    [~,Iend]=min(abs(xbins-Rmax));
else
    Iend=length(xbins);
end

SQ=zeros(1,length(Q));
if use_window==1
    for i=1:length(Q)
        SQ(i)=1+4*pi*pho*trapz(xbins(1:Iend),(pdf(1:Iend)-1).*sin(Q(i)*xbins(1:Iend)).*xbins(1:Iend).^2./(Q(i)*xbins(1:Iend)).*sin(pi*xbins(1:Iend)/Rmax)./(pi*xbins(1:Iend)/Rmax));
        
        %             SQal(i)=4*pi*concentrn*sqrt(ca*cb)*trapz(xbins(1:Iend),(pdf(1:Iend)-1).*sin(Q(i)*xbins(1:Iend)).*xbins(1:Iend).^2./(Q(i)*xbins(1:Iend)).*sin(pi*xbins(1:Iend)/Rmax)./(pi*xbins(1:Iend)/Rmax));
        %             SQfz(i)=SQal(i)/sqrt(ca*cb)+1;
        
    end
else
    for i=1:length(Q)
        SQ(i)=1+4*pi*pho*trapz(xbins(1:Iend),(pdf(1:Iend)-1).*sin(Q(i)*xbins(1:Iend)).*xbins(1:Iend).^2./(Q(i)*xbins(1:Iend)));
    end
end