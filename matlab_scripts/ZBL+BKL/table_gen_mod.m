clear, clc, close all;

%% Just for oxygen search for this line
%        pp=polyfit(xx(120:300),yy(120:300),5); 
%        figure;
%        plot(xx,yy,xx,polyval(pp,xx))
%        xlim([xx(120) xx(300)])
% for other interactions :
%        pp=polyfit(xx(90:300),yy(90:300),6); 
%        figure;
%        plot(xx,yy,xx,polyval(pp,xx))
%        xlim([xx(90) xx(300)])

%%
factor2=1/(4.184*1000/6.022140857e23/1.60217662e-19); % from eV to kcal/mol

cutoff=7.5;
table_cutoff=10;

fprintf('Mountjoy Potential\n');
%A= [101093.347  281375.657  13702.905*factor2  1844.7458*factor2];
%rho=[0.243838   0.195628    0.193817      0.343645];
%C  =[707.96963  737.8796173 54.681*factor2  192.58*factor2];

A=[ 639279.02  42523.2854]; %.*factor2; %2029.22, 1844.7458(old)
rho=[0.1819   0.3436];
C=[2003.02  4440.962203]; %.*factor2;
Atom1={ 'P'   'O' };
Atom2={ 'O'   'O' };
r_zbl=[   0.75   1     ]; %Pure ZBL before this
r_vdw=[   0.98   1.435 ]; %Pure buckingham after this
cut_zbl=[   20   20    ];
%r_zbl=[ 0.75  0.75   0.75   0.929  0.75   0.75  0.75  0.75];
%r_vdw=[ 1.2   1.12   1.1    1.495  1.35   1.4   1.1   1.1];
%cut_zbl=[20   20     20     20    20      20      20      20];

Z=    [ 8   15 ];
Aname={'O' 'P'};

for i=1:length(Atom1)
    if length(Atom1{i})~=2
        if length(Atom1{i})==1
            Atom1{i}=[' ',Atom1{i}];
        else
            error('check atom name in Atom1!');
        end
    end
    if length(Atom2{i})~=2
        if length(Atom2{i})==1
            Atom2{i}=[' ',Atom2{i}];
        else
            error('check atom name in Atom2!');
        end
    end
end
for i=1:length(Aname)
    if length(Aname{i})~=2
        if length(Aname{i})==1
            Aname{i}=[' ',Aname{i}];
        else
            error('check atom name in Aname!');
        end
    end
end

paircheck=cell2mat([Atom1',Atom2']);


npair=factorial(length(Z))/factorial(length(Z)-2)/factorial(2)+length(Z);
pair=cell(npair,2);
pZBL=zeros(npair,2);
pbuck=zeros(npair,3);
pcat=zeros(npair,3);
indbuck=1:length(Atom1);
l=1;
for i=1:length(Z)
    for j=i:length(Z)
        pair(l,:)={Aname{i}, Aname{j}};
        pZBL(l,:)=[Z(i), Z(j)];
        
        yesbuck=0;
        for k=1:size(paircheck)
            if strcmp(paircheck(k,:),[Aname{i}, Aname{j}]) || strcmp(paircheck(k,:),[Aname{j}, Aname{i}])
                yesbuck=k;
                break;
            end
        end
        if yesbuck>0
            pbuck(l,:)=[A(yesbuck) rho(yesbuck) C(yesbuck)];
            pcat(l,:)=[r_zbl(yesbuck) r_vdw(yesbuck) cut_zbl(yesbuck)];
        end
        fprintf('%s-%s: A%15.4f, rho%8.4f, C%10.4f\n',pair{l,1},pair{l,2},pbuck(l,:));
        l=l+1;
    end
end

for i=1:npair
    
    if pbuck(i,1)==0
    i   
        
        xxx=linspace(0.01,table_cutoff,10000);
        Ezbl=zbl(xxx,pZBL(i,1),pZBL(i,2));%,0.001,2);
        Ezbl(xxx>cutoff)=0;
        Ezbl(Ezbl<0)=0;
        
%         pp=spline(xxx,Ezbl);
%         Fpp=fnder(pp,1);
%         Fzbl=-ppval(Fpp,xxx');
%         Fzbl(xxx>cutoff)=0;
%         Fzbl(Ezbl==0)=0;
%         
                figure;
                plot(xxx,Ezbl,xxx,zbl(xxx,pZBL(i,1),pZBL(i,2)))
                title(['Energy, ',pair{i,1},pair{i,2}]);
                xlim([0.2 1])
        
                figure;
                plot(xxx,Fzbl,xxx(2:end)-(xxx(2)-xxx(1))/2,-diff(ppval(pp,xxx))./diff(xxx))
                title(['Force, ',pair{i,1},pair{i,2}]);
        
        pairname=[pair{i,1},pair{i,2}];
        pairname=pairname(pairname~=' ');
        fd=fopen(['BKS_ZBL_',pairname],'w');
        id=1:length(xxx);
        fprintf(fd,[pairname,'\nN %4.0f R %6.2f%6.2f\n\n'],[length(id) xxx(1) xxx(end)]);
        fprintf(fd,'%5.0f%25.15f%25.15e%25.15e\n',[id',xxx',Ezbl',Fzbl]');
        fclose(fd);
    elseif pbuck(i,3)==0
        xxx=linspace(0.01,table_cutoff,10000);
        Ezbl=buck(xxx,pbuck(i,1),pbuck(i,2),pbuck(i,3));
        Ezbl(xxx>cutoff)=0;
        Ezbl(Ezbl<0)=0;
        
        pp=spline(xxx,Ezbl);
        Fpp=fnder(pp,1);
        Fzbl=-ppval(Fpp,xxx');
        Fzbl(xxx>cutoff)=0;
        Fzbl(Ezbl==0)=0;
        
                figure;
                plot(xxx,Ezbl,xxx,buck(xxx,pbuck(i,1),pbuck(i,2),pbuck(i,3)))
                title(['Energy, ',pair{i,1},pair{i,2}]);
                xlim([0 1])
        
                figure;
                plot(xxx,Fzbl,xxx(2:end)-(xxx(2)-xxx(1))/2,-diff(ppval(pp,xxx))./diff(xxx))
                title(['Force, ',pair{i,1},pair{i,2}]);
        
        pairname=[pair{i,1},pair{i,2}];
        pairname=pairname(pairname~=' ');
        fd=fopen(['BKS_ZBL_',pairname],'w');
        id=1:length(xxx);
        fprintf(fd,[pairname,'\nN %4.0f R %6.2f%6.2f\n\n'],[length(id) xxx(1) xxx(end)]);
        fprintf(fd,'%5.0f%25.15f%25.15e%25.15e\n',[id',xxx',Ezbl',Fzbl]');
        fclose(fd);
    else
        
        xx=[0.01:0.005:pcat(i,1) pcat(i,2):0.005:cutoff*1.5];
        yy=[zbl(0.01:0.005:pcat(i,1),pZBL(i,1),pZBL(i,2)),buck(pcat(i,2):0.005:cutoff*1.5,pbuck(i,1),pbuck(i,2),pbuck(i,3))];
        
        pp=polyfit(xx(110:300),yy(110:300),5); 
        figure;
        plot(xx,yy,xx,polyval(pp,xx))
        xlim([xx(120) xx(300)])
        
        xx=[0.01:0.001:pcat(i,1)*0.7 pcat(i,1)*0.8:0.001:pcat(i,2)*0.95 pcat(i,2)*1.1:0.001:cutoff]';
        yy=[zbl(0.01:0.001:pcat(i,1)*0.7,pZBL(i,1),pZBL(i,2)),polyval(pp,pcat(i,1)*0.8:0.001:pcat(i,2)*0.95),buck(pcat(i,2)*1.1:0.001:cutoff,pbuck(i,1),pbuck(i,2),pbuck(i,3))];
        yy=smooth(yy,'rloess');
        
        pp=csape(xx,yy);
        xxx=linspace(0.01,table_cutoff,10000);
        yyy=ppval(pp,xxx');
        yyy(xxx>cutoff)=0;
        
        figure;
        plot(xxx,zbl(xxx,pZBL(i,1),pZBL(i,2)),xxx,buck(xxx,pbuck(i,1),pbuck(i,2),pbuck(i,3)),xxx,yyy);
        xlim([pcat(i,1)*0.5 pcat(i,2)*2]);
        ylim([-2 50000]);
        legend('ZBL',[pair{i,1},pair{i,2}],'cat');
        
        %         ylim([-2 2000])
        %         xlim([pcat(i,1)*0.7 pcat(i,2)*2]);
        
        r_check=pcat(i,1)*0.5:0.01:pcat(i,2)*1.1;
        figure;
        plot(r_check(2:end),diff(ppval(pp,r_check)))
        title('1st der');
        legend('ZBL',[pair{i,1},pair{i,2}],'cat');
        figure;
        plot(r_check(3:end),diff(diff(ppval(pp,r_check))))
        title('2nd der');
        legend('ZBL',[pair{i,1},pair{i,2}],'cat');
        
        Fpp=fnder(pp,1);
        Fzbl=-ppval(Fpp,xxx');
        Fzbl(xxx>cutoff)=0;
        
        
        pairname=[pair{i,1},pair{i,2}];
        pairname=pairname(pairname~=' ');
        fd=fopen(['BKS_ZBL_',pairname],'w');
        id=1:length(xxx);
        fprintf(fd,[pairname,'\nN %4.0f R %6.2f%6.2f\n\n'],[length(id) xxx(1) xxx(end)]);
        fprintf(fd,'%5.0f%25.15f%25.15e%25.15e\n',[id',xxx',yyy,Fzbl]');
        
        fclose(fd);
        
    end
    
    system(['cat BKS_ZBL_',pairname,'>> BKS.',num2str(cutoff),'.',num2str(table_cutoff),'.table.new']);
end