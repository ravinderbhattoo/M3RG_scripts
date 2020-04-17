clear,clc;

% ===== read lammps trajectory: =====
ftrj=fopen('./md2.lammpstrj');
[timestep,box,id,label,coordinate]=read_lammpstrj(ftrj,1);
fclose(ftrj);
save ../coors timestep box id label coordinate;

load ../coors;

% ==== atomic definiations and parameters for pdf ====
alabel=[  1   2    3   4    5];
anames={'Na','B','Si','O','Ca'};
cutoff=10;
dr=0.01;

frames=[length(timestep)];
% frames=sort([20:5:502, 21:5:503, 22:5:504]);
% frames=[500 501 502];

% Neutron scattering length
element={'O',   'Si', 'Al',  'Ca', 'Na' ,'K' ,'Li', 'B'};
bb=   [5.803 4.1491 3.449 4.700  3.63   3.67 -1.90  6.65]; % NIST http://www.ncnr.nist.gov/resources/n-lengths/, B from ???


nr=length(0:dr:cutoff)-1;
aconc=zeros(1,length(alabel));
n_pair=length(alabel)*(length(alabel)+1)*0.5;
tot_pdf=zeros(length(frames),nr);
pdf_group=zeros(n_pair,nr,length(frames));
for nnn=1:length(frames)
    pair=zeros(n_pair,2);
    cij=zeros(n_pair,1);
    bij=zeros(n_pair,1);
    jjj=1;
    for i=1:length(alabel)
        for ii=i:length(alabel)
            pair(jjj,:)=[alabel(i),alabel(ii)];

            bij(jjj)=bb(strcmp(anames{i},element))*bb(strcmp(anames{ii},element));

            aconc1=sum(label(:,frames(nnn))==alabel(i))/size(label,1);
            aconc2=sum(label(:,frames(nnn))==alabel(ii))/size(label,1);
            cij(jjj)=aconc1*aconc2;

            coorsc=coordinate(label(:,frames(nnn))==alabel(i),:,frames(nnn));
            coorss=coordinate(label(:,frames(nnn))==alabel(ii),:,frames(nnn));
            [xbins,pdf]=PDFcal(box(:,:,frames(nnn)),coorsc,coorss,cutoff,dr);
            pdf_group(jjj,:,nnn)=pdf;

            jjj=jjj+1;
        end
    end

    bij(pair(:,1)==pair(:,2))=0.5*bij(pair(:,1)==pair(:,2));
    cb=sum(bij.*cij);
    bij=repmat(bij,1,nr);
    cij=repmat(cij,1,nr);
    tot_pdf(nnn,:)=sum(pdf_group(:,:,nnn).*bij.*cij)/cb;
    %     pdfnao(frames(nnn),:)=pdf_group(3,:);

end

save PDF frames xbins tot_pdf pdf_group alabel anames cutoff nr pair box;

color=lines(size(pdf_group,3));
figure;
hold on;
for i=1:10:size(pdf_group,3)
    plot(xbins,pdf_group(8,:,i),'Color',color(i,:));
end

load ../coors;
load PDF;

Rmax=16;
dq=0.01;
qrange=[pi/(2*Rmax) 25];

SQ_OO=zeros(size(pdf_group,3),length(floor(qrange(1)/dq):floor(qrange(2)/dq)));
SQ_BO=zeros(size(pdf_group,3),length(floor(qrange(1)/dq):floor(qrange(2)/dq)));
SQ_BB=zeros(size(pdf_group,3),length(floor(qrange(1)/dq):floor(qrange(2)/dq)));
SQ_SiO=zeros(size(pdf_group,3),length(floor(qrange(1)/dq):floor(qrange(2)/dq)));
SQ_SiSi=zeros(size(pdf_group,3),length(floor(qrange(1)/dq):floor(qrange(2)/dq)));
SQ_tot=zeros(size(tot_pdf,1),length(floor(qrange(1)/dq):floor(qrange(2)/dq)));
for i=1:size(pdf_group,3)
    pho=size(label,1)/(box(1,:,i)*cross(box(2,:,i),box(3,:,i))');

    [Q,SQ]=SQfz(xbins,pdf_group(13,:,i),pho,qrange,dq,1,Rmax);
    Q_OO=Q;
    SQ_OO(i,:)=SQ;

    [Q,SQ]=SQfz(xbins,pdf_group(8,:,i),pho,qrange,dq,1,Rmax);
    Q_BO=Q;
    SQ_BO(i,:)=SQ;
    
    [Q,SQ]=SQfz(xbins,pdf_group(6,:,i),pho,qrange,dq,1,Rmax);
    Q_BB=Q;
    SQ_BB(i,:)=SQ;
    
    [Q,SQ]=SQfz(xbins,pdf_group(11,:,i),pho,qrange,dq,1,Rmax);
    Q_SiO=Q;
    SQ_SiO(i,:)=SQ;

    [Q,SQ]=SQfz(xbins,pdf_group(10,:,i),pho,qrange,dq,1,Rmax);
    Q_SiSi=Q;
    SQ_SiSi(i,:)=SQ;
    
    [Q,SQ]=SQfz(xbins,tot_pdf(i,:),pho,qrange,dq,1,Rmax);
    Q_tot=Q;
    SQ_tot(i,:)=SQ;
end
save SQ Q SQ_tot SQ_OO SQ_BO SQ_BB SQ_SiO SQ_SiSi Rmax dq qrange frames;

color=lines(size(SQ_tot,1));
figure;
hold on;
for i=1:10:size(SQ_SiSi,1)
    plot(Q,SQ_tot(i,:),'Color',color(i,:));
end



% [~,fitparam,~,~]=lorentzfit(Q(fitrange),SQfz(fitrange));
% fitFSDP = fitparam(1)./((Q(rangeFSDP) - fitparam(2)).^2 + fitparam(3)) + fitparam(4);
% HM=fitparam(2)/2;
% Lorentz(HM,fitparam)




% fh=fopen('test.xyz','w');
% outlabel=cell(size(label(1)));
% for i=1:size(label,1)
%     outlabel{i}=anames{label(i)==alabel};
% end
% write_xyz(fh,size(label,1),coordinate(:,:,frames),outlabel)
% fclose(fh);