function [CN,Nlist,varargout]=CNcal(box,coorsr,coorss,cutoff)

% [Nlist,CN,bllist,balist]=CNcal(box,coorsr,coorss,cutoff)
% Calculate coordination number, by Bu
% CN: number of coorss around coorsr, within cutoff
%
% last update: 05/31/2015
% Calculate bond angle list (5/31/15)

extra_d=3;
na=floor((cutoff+extra_d)/box(1,1));
nb=floor((cutoff+extra_d)/box(2,2));
nc=floor((cutoff+extra_d)/box(3,3));

fracr=coorsr/box;
fracr(fracr>1)=fracr(fracr>1)-1;
fracr(fracr<0)=fracr(fracr<0)+1;
coorsr=fracr*box;

fracs=coorss/box;
fracs(fracs>1)=fracs(fracs>1)-1;
fracs(fracs<0)=fracs(fracs<0)+1;
coorss=fracs*box;
fracs=[fracs, (1:size(coorss,1))'];

frac_ext=fracs;

temp=frac_ext;
temp1=frac_ext;
temp2=frac_ext;
for i=1:na
    temp1(:,1)=temp(:,1)+i;
    temp2(:,1)=temp(:,1)-i;
    frac_ext=[frac_ext; temp1; temp2];
end

temp=frac_ext;
temp1=frac_ext;
temp2=frac_ext;
for i=1:nb
    temp1(:,2)=temp(:,2)+i;
    temp2(:,2)=temp(:,2)-i;
    frac_ext=[frac_ext; temp1; temp2];
end

temp=frac_ext;
temp1=frac_ext;
temp2=frac_ext;
for i=1:nc
    temp1(:,3)=temp(:,3)+i;
    temp2(:,3)=temp(:,3)-i;
    frac_ext=[frac_ext; temp1; temp2];
end

extr_x=(cutoff+extra_d)/box(1,1)-na;
extr_y=(cutoff+extra_d)/box(2,2)-nb;
extr_z=(cutoff+extra_d)/box(3,3)-nc;


temp1=frac_ext(frac_ext(:,1)<extr_x,:);
temp1(:,1)=temp1(:,1)+na+1;
temp2=frac_ext(frac_ext(:,1)>1-extr_x,:);
temp2(:,1)=temp2(:,1)-na-1;
frac_ext=[frac_ext; temp1; temp2];

temp1=frac_ext(frac_ext(:,2)<extr_y,:);
temp1(:,2)=temp1(:,2)+nb+1;
temp2=frac_ext(frac_ext(:,2)>1-extr_y,:);
temp2(:,2)=temp2(:,2)-na-1;
frac_ext=[frac_ext; temp1; temp2];

temp1=frac_ext(frac_ext(:,3)<extr_z,:);
temp1(:,3)=temp1(:,3)+nc+1;
temp2=frac_ext(frac_ext(:,3)>1-extr_z,:);
temp2(:,3)=temp2(:,3)-nc-1;
frac_ext=[frac_ext; temp1; temp2];

coors_ext=frac_ext(:,1:3)*box;
% figure;
% plot(frac_ext(:,1),frac_ext(:,2),'.',frac(:,1),frac(:,2),'.')
% figure;
% plot(frac_ext(:,2),frac_ext(:,3),'.',frac(:,2),frac(:,3),'.')

cutoff=cutoff+0.001;
Nlist=cell(size(coorsr,1),1);
CN=zeros(size(coorsr,1),1);
if nargout >= 3
    bllist=cell(size(coorsr,1),1);
end

if nargout == 4
    idext=1:size(coors_ext,1);
    balist=cell(size(coorsr,1),1);
end
for i=1:size(coorsr,1)
    ref=coorsr(i,:);
    Icand=coors_ext(:,1)>ref(:,1)-cutoff & coors_ext(:,1)<ref(:,1)+cutoff &...
        coors_ext(:,2)>ref(:,2)-cutoff & coors_ext(:,2)<ref(:,2)+cutoff &...
        coors_ext(:,3)>ref(:,3)-cutoff & coors_ext(:,3)<ref(:,3)+cutoff;
    idcand=frac_ext(Icand,4);
    dcand=coors_ext(Icand,:)-repmat(coorsr(i,:),sum(Icand),1);
    ddcand=sqrt(dcand(:,1).^2+dcand(:,2).^2+dcand(:,3).^2);
    
    Inbr=ddcand<=cutoff & ddcand>0.001;
    Nlist{i}=idcand(Inbr);
    CN(i)=length(Nlist{i});
    
    if nargout >= 3
        bllist{i}=ddcand(Inbr);
    end
    
    if nargout == 4
        idcandext=idext(Icand);
        temp0=coorsr(i,:);
        coors_nbr=coors_ext(idcandext(Inbr),:);
        l=1;
        for jj=1:size(coors_nbr,1)
            tempa=coors_nbr(jj,:);
            for kk=jj+1:size(coors_nbr,1)
                tempb=coors_nbr(kk,:);
                balist{i}(l,1)=rad2deg(acos((tempa-temp0)*(tempb-temp0)'/norm(tempa-temp0)/norm(tempb-temp0)));
                l=l+1;
            end
        end
    end
end

if nargout >= 3
    varargout{1}=bllist;
end

if nargout == 4
    varargout{2}=balist;
end