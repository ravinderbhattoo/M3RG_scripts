function [xbins,pdf,varargout]=PDFcal(box,coorsr,coorss,cutoff,binwidth)
extra_d=10;
na=floor((cutoff+extra_d)/box(1,1));
nb=floor((cutoff+extra_d)/box(2,2));
nc=floor((cutoff+extra_d)/box(3,3));


bins=0:binwidth:cutoff;
xbins=(bins(1:end-1)+bins(2:end))/2;
% ddd=(cutoff+3)/min([box(1,1),box(2,2),box(3,3)]);

fracr=[coorsr/box, (1:size(coorsr,1))'];
fracr(fracr(:,1)>1,1)=fracr(fracr(:,1)>1,1)-1;
fracr(fracr(:,2)>1,2)=fracr(fracr(:,2)>1,2)-1;
fracr(fracr(:,3)>1,3)=fracr(fracr(:,3)>1,3)-1;
fracr(fracr(:,1)<0,1)=fracr(fracr(:,1)<0,1)+1;
fracr(fracr(:,2)<0,2)=fracr(fracr(:,2)<0,2)+1;
fracr(fracr(:,3)<0,3)=fracr(fracr(:,3)<0,3)+1;
fracs=[coorss/box, (1:size(coorss,1))'];
fracs(fracs(:,1)>1,1)=fracs(fracs(:,1)>1,1)-1;
fracs(fracs(:,2)>1,2)=fracs(fracs(:,2)>1,2)-1;
fracs(fracs(:,3)>1,3)=fracs(fracs(:,3)>1,3)-1;
fracs(fracs(:,1)<0,1)=fracs(fracs(:,1)<0,1)+1;
fracs(fracs(:,2)<0,2)=fracs(fracs(:,2)<0,2)+1;
fracs(fracs(:,3)<0,3)=fracs(fracs(:,3)<0,3)+1;

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
% plot3(coorss(:,1),coorss(:,2),coorss(:,3),'or');
% pbaspect([1 1 1])
% figure;
% plot3(coors_ext(:,1),coors_ext(:,2),coors_ext(:,3),'.b',coorss(:,1),coorss(:,2),coorss(:,3),'or')
% pbaspect([1 1 1])

% Nlist=cell(size(fracr,1),1);
% CN=zeros(size(fracr,1),1);
allpdf=zeros(size(coorsr,1),length(bins)-1);

if nargout == 7
    id_nbr=zeros(size(coorsr,1),10);
    dis_nbr=zeros(size(coorsr,1),10);
    x_nbr=zeros(size(coorsr,1),10);
    y_nbr=zeros(size(coorsr,1),10);
    z_nbr=zeros(size(coorsr,1),10);
elseif nargout ~= 2
    error('Need 5 variables for neighbour info!');
end

idext=(1:size(coors_ext,1))';
cutoff=cutoff+0.0001;
for i=1:size(coorsr,1)
    ref=coorsr(i,:);
    Icand=coors_ext(:,1)>ref(:,1)-cutoff & coors_ext(:,1)<ref(:,1)+cutoff &...
        coors_ext(:,2)>ref(:,2)-cutoff & coors_ext(:,2)<ref(:,2)+cutoff &...
        coors_ext(:,3)>ref(:,3)-cutoff & coors_ext(:,3)<ref(:,3)+cutoff;   
    dcand=coors_ext(Icand,:)-repmat(coorsr(i,:),sum(Icand),1);
    ddcand=[sqrt(dcand(:,1).^2+dcand(:,2).^2+dcand(:,3).^2),frac_ext(Icand,4),idext(Icand)];
    ddcand(ddcand(:,1)<=0.001,:)=[];
    [~,III]=sort(ddcand(:,1));

    if nargout == 7
        id_nbr(i,:)=ddcand(III(1:10),2)';
        dis_nbr(i,:)=ddcand(III(1:10),1)';
        x_nbr(i,:)=coors_ext(idext(III(1:10)),1)';
        y_nbr(i,:)=coors_ext(idext(III(1:10)),2)';
        z_nbr(i,:)=coors_ext(idext(III(1:10)),3)';
    end
    
    
    noatoms=histc(ddcand(:,1),bins);
    allpdf(i,:)=noatoms(1:end-1)'./(diff(bins)*4*pi.*xbins.^2 ...
        *(size(fracs,1)/(cross(box(1,:),box(2,:))*box(3,:)')));
end
pdf=mean(allpdf);
if nargout == 7
        varargout{3}=id_nbr;
        varargout{4}=dis_nbr;
        varargout{5}=x_nbr;
        varargout{6}=y_nbr;
        varargout{7}=z_nbr;
 end
