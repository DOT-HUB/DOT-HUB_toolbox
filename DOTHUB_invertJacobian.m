function invJ = DOTHUB_invertJacobian(hMesh,hBasis,SD,J_basis_wav1,J_basis_wav2,headMesh,lambda,RecMethod)

% 

%PRELIMINARIES
if isstruct(lambda)
    if ~isfield(lambda,'method')
        warning('Inversion method not provide, lambda.method can be ''Tikhonov'' (-DEFAULT), ''Cov'', ''CortexConst'' or ''SpatialRef'', default method used')
        lambda.method = 'Tikhonov';
    end
    if ~isfield(lambda,'value')
        warning('Regularization parameter not provide, using default lambda.value = 0.1;')
        lambda.value = 0.1;
    end
elseif isempty(lambda)
    warning('Regularization parameter not provide, using default Tikhonov regularization with lambda = 0.1;')
    lambda.value = 0.1;
    lambda.method = 'Tikhonov';    
else
    value = lambda; clear lambda;
    if numel(value) == 1
        lambda.value = value;
        lambda.method = 'Tikhonov';
    elseif numel(value) == 2
        lambda.value = value;
        lambda.method = 'SpatialReg'; 
    else
        error('Maximum number of regularization parameters is two.')
    end
end

%Wavelengths
wavelengths = SD.Lambda;

% Combine into matrix (wavelength x chromphore), convert units from cm-1/M to mm-1/uM
E1 = GetExtinctions(wavelengths(1));
E2 = GetExtinctions(wavelengths(2));
Eall = [E1(1:2)./1e7; E2(1:2)./1e7]; %2x2

if strfind(RecMethod,'stnd')
     if strfind(lambda.method,'Tikhonov')
        % Lambda has only one value so Tikhonov regularization is assumed
        
        
        
        % Inversion matrix, wavelentgh 1
        JJT = J_basis_wav1*J_basis_wav1';
        S=svd(JJT);
        invJ1 = J_basis_wav1' / (JJT + eye(length(JJT))*(lambda.value*max(S)));
        
        JJT = J_basis_wav2*J_basis_wav2';
        S=svd(JJT);
        invJ2 = J_basis_wav2' / (JJT + eye(length(JJT))*(lambda.value*max(S)));
        
        a = Eall(1,1);
        b = Eall(1,2);
        c = Eall(2,1);
        d = Eall(2,2);
        D = 1/(a*d - b*c);
%         if abs(D)<1e8
            invJ = D*[ d*invJ1 -b*invJ2;
                      -c*invJ1  a*invJ2];
%         else
%             % NOTE: ideally, wavelengths should have been selected better.
%             % Alternatively, try to regularize the matrix (in progress...)
%             error('Ill-condition matrix, aborting reconstruction')
%         end
        
     end   
    
elseif strfind(RecMethod,'multispec')
    %Building extinction coefficient dummy
    Ei = ones(size(J_basis_wav1));
    E = [];
    for c = 1:2 %Chromophore
        El = [];
        for i = 1:length(wavelengths)
            El = [El; Ei.*Eall(i,c)];
        end
        E = [E El];
    end

    %BUILD MULTISPECTRAL JACOBIAN
    fprintf('Building Multispectral Jacobian\n');
    Ji = [J_basis_wav1 J_basis_wav1; J_basis_wav2 J_basis_wav2]; %Tile wavelength-specific jacobians
    Jmulti = Ji.*E; % This has units of (d(ln(data/reference))/d(mm-1))*(mm-1/micromolar)


    disp('Inverting J ##############################################');

    if strfind(lambda.method,'Tikhonov')
        % Lambda has only one value so Tikhonov regularization is assumed

        % Reconstruction
        JJT = Jmulti*Jmulti';
        S=svd(JJT);
        invJ = Jmulti' / (JJT + eye(length(JJT))*(lambda.value*max(S)));

    elseif strfind(lambda.method,'Cov')

        nNodesBasis = length(J_basis_wav1);

        sigma_v = cov(ref_OD);
        sigma_u = sparse(1:2*nNodesBasis,1:2*nNodesBasis,1);
        JJT = Jmulti*sigma_u*Jmulti';
        lambda1 = lambda.value*trace(JJT)/trace(sigma_v);

        invJ = sigma_u*Jmulti' / ( JJT + sigma_v*lambda1);

    elseif strfind(lambda.method,'CortexConst')

        Jmulti_backup = Jmulti;
        Mask = zeros(size(headMesh.node,1),1);
        if max(headMesh.node(:,4)==5)
             %Assumes (1) Scalp, (2) Skull, (3) CSF, (3) GM, (5) WM
            Mask(headMesh.node(:,4)>3)=1;
        elseif max(headMesh.node(:,4)==4)
             %Assumes (1) ECT, (2) CSF, (3) GM, (4) WM
            Mask(headMesh.node(:,4)>2)=1;
        else
            error('Nodal tissue indices are not as expect for infant (1:4) or adult (1:5)');
        end
        bMask = hBasis.Map('M->S',Mask);
        bMask_nobrain = find(bMask==0);
        opc = 5;
        switch opc
            case 1        % Option 1: set to absolute minimum indicating change is minimal
                Jmulti = Jmulti_backup;

                tmp  = Jmulti(:,1:end/2);
                [val,ind] = min(abs(tmp(:)));
                Jmulti(:,bMask_nobrain) = sign(tmp(ind))*val;

                tmp  = Jmulti(:,end/2+1:end);
                [val,ind] = min(abs(tmp(:)));
                Jmulti(:,bMask_nobrain + hBasis.slen()) = sign(tmp(ind))*val;

            case 2        % Option 2: Deleting rows belonging to no-brain regions
                Jmulti = Jmulti_backup;

                Jmulti(:,bMask_nobrain + hBasis.slen()) = [];
                Jmulti(:,bMask_nobrain) = [];

            case 3       % Option 3: Setting no-brain related colums to zero
                Jmulti = Jmulti_backup;

                Jmulti(:,bMask_nobrain) = 0;
                Jmulti(:,bMask_nobrain + hBasis.slen()) = 0;

            case 4        % Option 4: Setting no-brain related colums to the absolute minimum of the column-mean
                m = mean(Jmulti,1);
                Jmulti = Jmulti_backup;

                tmp = m(1:end/2);
                [val,ind] = min(abs(tmp));
                Jmulti(:,bMask_nobrain) = sign(tmp(ind))*val;

                tmp = m(end/2+1:end);
                [val,ind] = min(abs(tmp));
                Jmulti(:,bMask_nobrain + hBasis.slen()) = sign(tmp(ind))*val;

            case 5        % Option 5: Setting no-brain related colums to the mean of the column-mean
                Jmulti = Jmulti_backup;

                tmp  = Jmulti(:,1:end/2);
                Jmulti(:,bMask_nobrain) = mean(tmp(:));

                tmp  = Jmulti(:,end/2+1:end);
                Jmulti(:,bMask_nobrain + hBasis.slen()) = mean(tmp(:));
        end

        JJT = Jmulti*Jmulti';
        S=svd(JJT);
        invJ = Jmulti' / (JJT + eye(length(JJT))*(lambda.value*max(S)));

    elseif strfind(lambda.method,'SpatialReg')
        % If lambda has two elements, it assumes spatial regularization as in:
        % White, B. (2012). Developing High-Density Diffuse Optical Tomography
        % for Neuroimaging. Washington University. pg 23, PhD Thesis.
        % (https://openscholarship.wustl.edu/etd/665/) Accesed on 16.04.2019

        lambda1 = lambda.value(1); % Typical Tikhonov regularization parameter
        lambda2 = lambda.value(2); % Spatial regularization parameter

        % Reconstruction
        JJT = Jmulti*Jmulti';  % Prepare Jacobian matrix for inversion (i.e. create square matrix)
        L = sqrt(diag(JJT) + lambda2*max(diag(JJT))); % Apply regularization
        Linv = 1./L; % Invert matrix
        %% find Atild
        Atild = zeros(size(Jmulti));
        for ind = 1:length(Linv)
            Atild(ind,:) = Jmulti(ind,:)*Linv(ind);
        end
        atildtatild = Atild*Atild';
        [satild] = svd(atildtatild);
        mxsatild = max(satild);

        % Apply spatial regularization
          val2binv = atildtatild;
          for ind = 1:length(satild)
              val2binv(ind,ind) = atildtatild(ind,ind) + lambda1*mxsatild;
          end

    %     %%Checks (verify regularization threshols)
    %     [satild_reg] = svd(val2binv);
    %     S=svd(JJT);
    %     figure;semilogy(satild,'linewidth',2),hold on,grid on,xlabel('Magnitude of singular values'),ylabel('Singular value index')
    %     Xlim = get(gca,'Xlim');plot(Xlim,lambda1*mxsatild*[1 1],'r','linewidth',2),semilogy(satild_reg,'k','linewidth',2)
    %     plot(S,'b--','linewidth',2),plot(Xlim,lambda.value(1)*max(S)*[1 1],'r--','linewidth',2),plot(S+lambda.value(1)*max(S),'k--','linewidth',2)
    %     hl= legend('Singular values (Spatial Reg)',['Regularization value \lambda\timesmax(\sigma) = ',num2str(lambda1*mxsatild,3),' (Spatial Reg)'],'Regularized Jacobian (Spatial Reg)',...
    %            'Singular values (Tikhonov Reg)',['Regularization value \lambda\timesmax(\sigma) = ',num2str(lambda.value(1)*max(S),3), '(Tikhonov Reg)'],'Regularized Jacobian (Tikhonov Reg)');
    %     set(hl,'Position',[0.1330    0.1109    0.6346    0.2440]),legend('boxoff')
    %     title(['\lambda_1 = ',num2str(lambda1),', \lambda_2 = ', num2str(lambda2)])

        clear JJT;   % Clear huge matrices for efficiency

        % Invert matrix
          inva = val2binv\Atild;
          clear Atild val2binv  % Clear huge matrices for efficiency
          invJ = zeros(size(inva))';
          for ind = 1:size(inva,1)
              invJ(:,ind) = inva(ind,:)*Linv(ind);
          end
    else
        error('Unknown method')
    end
end