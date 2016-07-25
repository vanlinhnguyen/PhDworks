close all; clc;
addpath('../../missingfuncs/isotropic');

%% INITIAL PARAMETERS
Nh=384; % Number of grid-point in HR field
k_max=192;

%% Loop all file
Nl=96; % Number of grid-point in HR field
nc1=netcdf('/data/PhDworks/isotropic/refdata_downsampled4.nc','clobber');
nc1('Nt')=0; 
nc1('Nx')=Nl;
nc1('Ny')=Nl;
nc1('Nz')=Nl;
nc1{'velocity_x'}=ncfloat('Nt','Nx','Ny','Nz');
nc1{'velocity_y'}=ncfloat('Nt','Nx','Ny','Nz');
nc1{'velocity_z'}=ncfloat('Nt','Nx','Ny','Nz');
count=0;

for file_ids=20:56
    start=tic();
    count=count+1;
    
    filename=strcat('FIELD-',num2str(file_ids,'%.3d'),'.nc');
    cmdString = ['scp nguyen@lmlm6-3.univ-lille1.fr:/data4/data/',filename,' /data/ISOTROPIC/data/'];
    fprintf(strcat('Loading:',filename,'\n'));
    [status, ~] = unix(cmdString);
    if status==0
        fprintf([filename,' downloaded in ',num2str(toc(start),'%.2f'),' seconds! \n']);

        nc2=netcdf(['/data/ISOTROPIC/data/', filename],'r'); 
        velocity_x_resolved=nc2{'velocity_x'}(:,:,:);
        velocity_y_resolved=nc2{'velocity_y'}(:,:,:);
        velocity_z_resolved=nc2{'velocity_z'}(:,:,:);
        close(nc2)

        cmdString = ['rm ',strcat('/data/ISOTROPIC/data/',filename)];
        fprintf(strcat('Delete:',filename,'\n'));
        [status, ~] = unix(cmdString);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% RATIO 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        start=tic();
        
        ratio=4;
        k_cutoff=k_max/ratio;
        Nl=2*k_cutoff; % Number of grid-point in LR field

        %% U component
        F = fftn(velocity_x_resolved);

        % Sample: remove high frequency
        ids=1:2*k_max;
        ids([k_cutoff+1:k_max k_max+1:2*k_max-k_cutoff])=[];
        F=F(ids,ids,ids);

        % Remove 3D wave number > k_cutoff
        k_1D_LR=[0:k_cutoff -k_cutoff+1:1:-1];
        
        for m=1:2*k_cutoff
            for n=1:2*k_cutoff
                for p=1:2*k_cutoff
                    if round(sqrt(k_1D_LR(m)^2+k_1D_LR(n)^2+k_1D_LR(p)^2))>k_cutoff
                        F(m,n,p)=0;
                    end
                end
            end
        end

        % ifftn back to physical space (normalize by size)
        velocity_x_downsampled=(Nl/Nh)^3*ifftn(F);
        clearvars F
        
        %% V component
        F = fftn(velocity_y_resolved);

        % remove high frequency
        ids=1:2*k_max;
        ids([k_cutoff+1:k_max k_max+1:2*k_max-k_cutoff])=[];
        F=F(ids,ids,ids);

        % Remove 3D wave number > k_cutoff
        k_1D_LR=[0:k_cutoff -k_cutoff+1:1:-1];
        for m=1:2*k_cutoff
            for n=1:2*k_cutoff
                for p=1:2*k_cutoff
                    if round(sqrt(k_1D_LR(m)^2+k_1D_LR(n)^2+k_1D_LR(p)^2))>k_cutoff
                        F(m,n,p)=0;
                    end
                end
            end
        end

        % ifftn back to physical space (normalize by size)
        velocity_y_downsampled=(Nl/Nh)^3*ifftn(F);
        clearvars F    
        
        
        %% W component
        F = fftn(velocity_z_resolved);

        % remove high frequency
        ids=1:2*k_max;
        ids([k_cutoff+1:k_max k_max+1:2*k_max-k_cutoff])=[];
        F=F(ids,ids,ids);

        % Remove 3D wave number > k_cutoff
        k_1D_LR=[0:k_cutoff -k_cutoff+1:1:-1];
        for m=1:2*k_cutoff
            for n=1:2*k_cutoff
                for p=1:2*k_cutoff
                    if round(sqrt(k_1D_LR(m)^2+k_1D_LR(n)^2+k_1D_LR(p)^2))>k_cutoff
                        F(m,n,p)=0;
                    end
                end
            end
        end

        % ifftn back to physical space (normalize by size)
        velocity_z_downsampled=(Nl/Nh)^3*ifftn(F);
        clearvars F   
        
        %% SAVE TO FILE        
        nc1{'velocity_x'}(count,:,:,:) = real(velocity_x_downsampled); 
        nc1{'velocity_y'}(count,:,:,:) = real(velocity_y_downsampled); 
        nc1{'velocity_z'}(count,:,:,:) = real(velocity_z_downsampled); 
        
        fprintf(strcat('Downsample by 4 in:',num2str(toc(start),'%.2f'),' seconds! \n'));

        
        
        %% Clear vars
        clearvars velocity_x_resolved velocity_y_resolved velocity_z_resolved

    else
        fprintf(strcat(filename,' not exist! \n'));
    end
end
close(nc1);