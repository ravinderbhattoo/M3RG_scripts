clear, clc;

% 11/12/2015, Bu, modify input method

format='dat';  % 'dat' for LAMMPS, 'gro' for gromacs, 'config' for DL_POLY

N=500;  % Total number of atoms

info_atoms=[1 22.990   0.6 % Na %format is id atomic_mass charge
    2 30.974  3.000   % P
    3 16.000  -1.2]; % O

% mole of oxide
% 0 for oxygen
molOXIDES=[1
    1
    0];

% Oxide*fOX2AT=ATOMS -> Number of cations in each oxide, 
% for example: 2 in Al2O3
% 0 for oxygen
fOX2AT=[2
    2
    0];

% mole% of cations, not oxide; e.g., 0.33 Na2O -> 0.33*2 Na
% 0 for oxygen and unavailable cations
molSPECIES = molOXIDES.*fOX2AT;


idO=info_atoms(:,3)<0;
if sum(idO)~=1
    error('Multiple anions!')
end
qO     =info_atoms(idO,3);
VALENCE=info_atoms(~idO,3)';
molCAT=molSPECIES(~idO)';


nCations=2*round(N*molCAT./sum(molCAT.*(1+VALENCE/abs(qO)))/2);
if sum(nCations(molCAT==0))~=0
    error('Check nCation calculation!');
end
nO=round(sum(nCations.*VALENCE)/abs(qO));
newN=nO+sum(nCations);
fprintf('Total number of atom:%5u (Targeted:%5u)\n',newN,N)
printstr=repmat(['%6.2f'],1,length(molCAT));
fprintf('                ')
fprintf('    %2.0f',info_atoms(~idO,1))
fprintf(['\nmole%% of cations:',printstr,'\n'],100*nCations/sum(nCations))
fprintf('                ')
fprintf('    %2.0f',info_atoms(~idO,1))
fprintf(['\nmole%% of oxides :',printstr,'\n'],100*nCations./fOX2AT(~idO)'/sum(nCations./fOX2AT(~idO)'))

a=nO^(1/3)*2*1.6; % estimate the size of the box
nx=ceil(nO^(1/3));
dx=a/nx;

coors_O=zeros(nx*nx*nx,3);
coors_O(:,1)=reshape(ones(1,nx*nx)'*(0:dx:(nx-1)*dx),[],1);
coors_O(:,2)=reshape(reshape(ones(1,nx)'*(0:dx:(nx-1)*dx),[],1)*ones(1,nx),[],1);
coors_O(:,3)=reshape((0:dx:(nx-1)*dx)'*ones(1,nx*nx),[],1);

coors_cat=coors_O+ones(size(coors_O,1),1)*[dx/2 dx/2 dx/2];

rng('shuffle');
coors_O=coors_O(randperm(size(coors_O,1),nO),:);
coors_cat=coors_cat(randperm(size(coors_cat,1),sum(nCations)),:);


%%%%% LAMMPS DAT %%%%%
if strcmp(format,'dat')
    all_CAT_ID=zeros(sum(nCations),1);
    kk=0;
    for i=1:size(info_atoms,1)-1
        all_CAT_ID(kk+1:kk+nCations(i))=repmat(sortrows(info_atoms(i,1),1),nCations(i),1);
        kk=kk+nCations(i);
    end
    fdat=fopen('glass.dat','w');
    write_dat(fdat,'charge',sortrows(info_atoms,1),[a 0 0; 0 a 0; 0 0 a],[coors_cat;coors_O],[all_CAT_ID;repmat(info_atoms(idO,1),size(coors_O,1),1)]);
    
end

% %%%%% GROMACS GRO %%%%%
% if strcmp(format,'gro')
    
% CAT_names={'   '};
% for i=1:length(CATIONS)
%     CAT_names=strcat(CAT_names,{'   '},CATIONS{i});
% end
% CATIONS={'Na','B','Si'};
% RES    ={  'NAN', 'SIS'};  % for GROMACS
% if strcmp(format,'gro')
% CAT_label=cell(sum(nCations),1);
% RES_label=CAT_label;
% kk=0;
% for i=1:length(CATIONS)
%     CAT_label(kk+1:kk+nCations(i))=mat2cell(repmat(CATIONS{i},nCations(i),1),ones(nCations(i),1),2);
%     RES_label(kk+1:kk+nCations(i))=mat2cell(repmat(RES{i},nCations(i),1),ones(nCations(i),1),3);
%     kk=kk+nCations(i);
% end


%     gro=fopen('glass.gro','w');
%     
%     fprintf(gro,'%s   Glass structure for melting, %7u atoms\n',CAT_names{1},N);
%     fprintf(gro,'%5u\n',N);
%     
%     for i=1:sum(nCations)
%         fprintf(gro,'    1%3s     %2s%5u%8.3f%8.3f%8.3f\n',RES_label{i},CAT_label{i},i,coors_cat(i,:)*0.1);
%     end
%     for i=1:nO
%         fprintf(gro,'    1OXY     O %5u%8.3f%8.3f%8.3f\n',i+sum(nCations),coors_O(i,:)*0.1);
%     end
%     
%     fprintf(gro,'%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f',0.1*[a a a 0 0 0 0 0 0]');
% end
% 
% end
%
% %%%%%  DL_POLY CONFIG %%%%%
% if strcmp(format,'config')
%     cfg=fopen('CONFIG','w');
%     fprintf(cfg,'%s   Glass structure for melting, %7u atoms\n',CAT_names{1},N);
%     fprintf(cfg,'   0    1\n');
%     fprintf(cfg,'%15.10f%15.10f%15.10f\n%15.10f%15.10f%15.10f\n%15.10f%15.10f%15.10f\n',[a 0 0; 0 a 0;0 0 a]');
%     
%     for i=1:sum(nCations)
%         fprintf(cfg,'%2s      %5u\n%15.10f%15.10f%15.10f\n',CAT_label{i},i,coors_cat(i,:));
%     end
%     
%     
%     for i=1:nO
%         fprintf(cfg,' O      %5u\n%15.10f%15.10f%15.10f\n',i+sum(nCations),coors_O(i,:));
%     end
% end

fclose all;