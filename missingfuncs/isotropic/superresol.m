%SUPERRESOL super-resolution using coupled dictionaries (SR1, SR2 or SR3)
%
%IN:
%   X_LR - input LR field
%   {D_HR,D_LR} - coupled dictionaries (trained a-priori)
%   mode - to choose SR1, SR2 or SR3 accordingly
%   params - parameters to do SR (see SPAMS documentation)
%   varargin - optional arguments for when mode = 'SR3' (derivetive kernels
%   and V_pca for dimension reduction)

%OUT:
%   X_HR - output HR field
%   CoefMatrix  - coefficient of sparse representation

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function [X_HR,CoefMatrix] = superresol(X_LR, D_HR, D_LR, params, mode, grids, sizes, varargin)
    switch(mode)
        case('SR1')
            patchsize_h = sqrt(size(D_HR, 1)); 
            patchsize_l = sqrt(size(D_LR, 1));
            gridz_l=grids{1}; gridy_l=grids{2}; gridz_h=grids{3}; gridy_h=grids{4};        
            Nz_h=sizes{3}; Ny_h=sizes{4}; % Nz_l=sizes{1}; Ny_l=sizes{2};
            
            % collect patches
            patches_LR=extract_patches(X_LR, {gridz_l,gridy_l}, patchsize_l);
            m_LR=mean(patches_LR,1); norm_LR=sqrt(sum(patches_LR.^2,1));
            patches_LR=patches_LR - repmat(m_LR,patchsize_l^2,1);
            patches_LR=patches_LR./repmat(norm_LR,patchsize_l^2,1);
            
            % solve sparse code and reconstruct HR patches 
            CoefMatrix=mexLasso(patches_LR,D_LR,params);
            rec_LR=D_LR*CoefMatrix; 
%             fprintf(['NRMSE of feature reconstruction:',num2str(sqrt(sum((patches_LR(:)-rec_LR(:)).^2))/sqrt(sum(patches_LR(:).^2)),'%.3f'),'\n']);

            patches_HR = D_HR*CoefMatrix;
            patches_HR = patches_HR.*repmat(norm_LR,patchsize_h^2,1);
            patches_HR = patches_HR + repmat(m_LR, patchsize_h^2, 1);

            X_HR=putback_patches(patches_HR, {gridz_h,gridy_h}, {Nz_h, Ny_h});
        
        case('SR2')
            patchsize_h = sqrt(size(D_HR, 1)); 
            gridz_h=grids{3}; gridy_h=grids{4};        
            Nz_h=sizes{3}; Ny_h=sizes{4};
            
            % collect patches
            patches_HF=extract_patches(X_LR, {gridz_h,gridy_h}, patchsize_h);
            m_HF=mean(patches_HF,1); norm_HF=sqrt(sum(patches_HF.^2,1));
            patches_HF=patches_HF - repmat(m_HF,patchsize_h^2,1);
            patches_HF=patches_HF./repmat(norm_HF,patchsize_h^2,1);
            
            % solve sparse code and reconstruct HR patches 
            CoefMatrix=mexLasso(patches_HF,D_LR,params);
            rec_HF=D_LR*CoefMatrix; 
%             fprintf(['NRMSE of feature reconstruction:',num2str(sqrt(sum((patches_HF(:)-rec_HF(:)).^2))/sqrt(sum(patches_HF(:).^2)),'%.3f'),'\n']);
%             fprintf(['Average number of nonzero coefficients ',num2str(round(nnz(CoefMatrix)/size(CoefMatrix,2)),'%2d'),'\n']);

            patches_LF = D_HR*CoefMatrix;
            patches_LF = patches_LF.*repmat(norm_HF,patchsize_h^2,1);
            patches_LF = patches_LF + repmat(m_HF, patchsize_h^2, 1);

            X_HR=putback_patches(patches_LF, {gridz_h,gridy_h}, {Nz_h, Ny_h});
            
            
        case('SR3')
            if size(varargin,2) == 0
                error('Error: No gradient filter found!')
            end
            
            if length(grids)~=2
                error('Pass 2 grids at HR only')    
            else
                gridz=grids{1};
                gridy=grids{2};
            end
            
            filters=varargin{1};
            V_pca=varargin{2};
            
            patchsize = sqrt(size(D_HR, 1)); % same for HR and LR
            [Ny, Nz] = size(X_LR);    
                        
            noFil=length(filters);
            
            % compute gradient
            X_LR_Fea=zeros(Ny,Nz,noFil);
            for kk=1:noFil
                X_LR_Fea(:, :, kk) = conv2(X_LR, filters{kk}, 'same');
            end 
                
            % prepare to reconstruct
            X_HR = zeros(Ny, Nz); 
            cntMat = zeros(Ny, Nz); 

            % collect patches
            patches_LR_Fea=zeros(noFil*patchsize*patchsize,length(gridz)*length(gridy)); 
            for ii = 1:length(gridz)
                for jj = 1:length(gridy)
                    zz = gridz(ii);
                    yy = gridy(jj);
                     
                    patch_LR_Fea = zeros(noFil*patchsize*patchsize,1);
                    for kk=1:noFil
                        patch_LR_Fea((kk-1)*patchsize*patchsize+1:kk*patchsize*patchsize) = X_LR_Fea(yy:yy+patchsize-1, zz:zz+patchsize-1,kk);
                    end
                     
                    patches_LR_Fea(:,(ii-1)*length(gridy)+jj)=patch_LR_Fea;
                end
            end
            patches_LR_Fea_pca = V_pca' * patches_LR_Fea;
            
            m_patchesLR_Fea_pca=mean(patches_LR_Fea_pca,2);
            for i=1:size(patches_LR_Fea_pca,2)
                patches_LR_Fea_pca(:,i) = patches_LR_Fea_pca(:,i) - m_patchesLR_Fea_pca;
            end
            norm_patchesLR_Fea_pca=sqrt(sum(patches_LR_Fea_pca.^2));
            for i=1:size(patches_LR_Fea_pca,1)
                patches_LR_Fea_pca(i,:) = patches_LR_Fea_pca(i,:)./norm_patchesLR_Fea_pca;
            end
            
            CoefMatrix=mexLasso(patches_LR_Fea_pca,D_LR,params);
            
            rec_Fea=D_LR*CoefMatrix; 
%             fprintf(['NRMSE of feature reconstruction:',num2str(sqrt(sum((patches_LR_Fea_pca(:)-rec_Fea(:)).^2))/sqrt(sum(patches_LR_Fea_pca(:).^2)),'%.3f'),'\n']);
            patches_HR = D_HR*CoefMatrix;


            % put HR patches back
            for ii = 1:length(gridz) 
                for jj = 1:length(gridy)
                    zz = gridz(ii);
                    yy = gridy(jj);
                    patch_HR = reshape(patches_HR(:,(ii-1)*length(gridy)+jj),[patchsize, patchsize]);                    
                    X_HR(yy:yy+patchsize-1, zz:zz+patchsize-1) = X_HR(yy:yy+patchsize-1, zz:zz+patchsize-1) + patch_HR;
                    cntMat(yy:yy+patchsize-1, zz:zz+patchsize-1) = cntMat(yy:yy+patchsize-1, zz:zz+patchsize-1) + 1;
                end
            end
            % fill in the empty with bicubic interpolation
            idx = (cntMat < 1); 
            cntMat(idx) = 1; X_HR(idx) = 0;
            X_HR = X_HR./cntMat;        
             
        otherwise
            fprintf('mode can be SR1, SR2 or SR3 only!\n');
    end
end
