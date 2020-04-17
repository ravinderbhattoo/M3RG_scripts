function [timestep,box,data]=read_dump_peri(fd,flag_sort,n_column)

% [timestep,box,id,label,coordinate]=read_lammpstrj(fileID,if_sort)
% Read LAMMPS trajectory file, by Bu and Alex
%
% last updated: 05/27/2015
% Automatically deal with ITEM line (5/27/15)
%
% The LAMMPS trajectory is indicated by the file handle, fileID. Use fopen
% to open the file and obtain the fileID value. When you finish reading,
% close the file by calling fclose(fileID).
%
% OUTPUT:
% 1. timestep, timestep of each frame. Array of length nframes, nframes is
% the number of frames.
% 2. box, box vectors of each frame. (3 x 3 x nframes) matrix.
% 3. id, indice of atoms for each frame. (natoms x nframes) matrix, natoms
% is the number of atoms.
% 5. label, labels of atoms for each frame. (natoms x nframes) matrix.
% 6. coordinates, coordiantes of atoms for each frame.
% (natoms x 3 X nframes) matrix.

fgetl(fd);
fgetl(fd);
fgetl(fd);
na=str2double(fgetl(fd));
temps=fgetl(fd);

if ~strcmp(temps,'ITEM: BOX BOUNDS xy xz yz pp pp pp') || ~strcmp(temps,'ITEM: BOX BOUNDS xy xz yz ss ss ss')
    fgetl(fd);
    fgetl(fd);
    fgetl(fd);
    ITEMline=fgetl(fd);
    frewind(fd);
    datastring=[repmat('%f',1,n_column),'\n'];
    %frame=textscan(fd,['ITEM: TIMESTEP\n%f\nITEM: NUMBER OF ATOMS\n%f\nITEM: BOX BOUNDS pp pp pp\n%f%f\n%f%f\n%f%f\n',ITEMline,'\n',repmat(datastring,1,na)]);
    frame=textscan(fd,['ITEM: TIMESTEP\n%f\nITEM: NUMBER OF ATOMS\n%f\nITEM: BOX BOUNDS ss ss ss\n%f%f\n%f%f\n%f%f\n',ITEMline,'\n',repmat(datastring,1,na)]);    
    nframe=length(frame{1});
    timestep=frame{1};
    
    temp=cell2mat(frame);
    data=permute(reshape(temp(:,9:end)',n_column,na,[]),[2,1,3]);
    
    box=zeros(3,3,nframe);
    box(1,1,:)=frame{4}-frame{3};
    box(2,2,:)=frame{6}-frame{5};
    box(3,3,:)=frame{8}-frame{7};
    
    if flag_sort==1
        for i=1:size(data,3)
            [~,I]=sort(data(:,1,i));
            data(:,:,i)=data(I,:,i);
        end
        if diff(data(:,1,i),1,2)~=0
            error('check sorting!')
        end
    end
else
    fgetl(fd);
    fgetl(fd);
    fgetl(fd);
    ITEMline=fgetl(fd);
    frewind(fd);
    
    datastring=[repmat('%f',1,n_column),'\n'];
    frame=textscan(fd,['ITEM: TIMESTEP\n%f\nITEM: NUMBER OF ATOMS\n%f\nITEM: BOX BOUNDS xy xz yz ss ss ss\n%f%f%f\n%f%f%f\n%f%f%f\n',ITEMline,'\n',repmat(datastring,1,na)]);
    
    nframe=length(frame{1});
    timestep=frame{1};
    
    data=permute(reshape(cell2mat(frame(12:end)),n_column,na,[]),[2,1,3]);
    
    box=zeros(3,3,nframe);
    xy=frame{5}';
    xz=frame{8}';
    yz=frame{11}';
    xlo=frame{3}-min([zeros(1,nframe);xy;xz;xy+xz])';
    xhi=frame{4}-max([zeros(1,nframe);xy;xz;xy+xz])';
    ylo=frame{6}-min([zeros(1,nframe);yz])';
    yhi=frame{7}-max([zeros(1,nframe);yz])';
    
    box(1,1,:)=xhi-xlo;
    box(2,1,:)=xy;
    box(2,2,:)=yhi-ylo;
    box(3,1,:)=xz;
    box(3,2,:)=yz;
    box(3,3,:)=frame{10}-frame{9};
    
    if flag_sort==1
        for i=1:size(data,3)
            [~,I]=sort(data(:,1,i));
            data(:,:,i)=data(I,:,i);
        end
        if diff(data(:,1,i),1,2)~=0
            error('check sorting!')
        end
    end
end


