function out = retina_V1_model(target, background, tx, ty, fx, fy, Parameters, Stacks)

%retina_V1_model predicts the contrast detection threshold of a target in a background. The location of the target is at
%(tx,ty) degrees from the center of the background, and the fixation point is at (fx,fy) degrees from the center of the
%background. The center of any n x m image is defined to be at (floor(n/2)+1,floor(m/2)+1). Pixels per degree is fixed 
%at 120. 

%REQUIRED INPUTS:
%target = grayscale target image, range: [-127,127]. Note that the target is defined as a contrast image that can be 
%         added to the background, which has a range of [0,255]. When added to the background, [-127,127] maps to 
%         [1,255]. Thus, zero maps to 128.
%background = grayscale background image, range: [0,255]. Each dimension of the background image must be at least twice
%             the largest dimension of the target image.
%tx = x coordinate (in degrees) of target image center relative to center of background.
%ty = y coordinate (in degrees) of target image center relative to center of background.
%fx = x coordinate (in degrees) of fixation location relative to center of background.
%fy = y coordinate (in degrees) of fixation location relative to center of background.

%OPTIONAL INPUTS:
%Parameters = a structure whose fields specify the parameters of the RV1 model. To use default parameter values, use [].
%             For example: X = retina_V1_model(target,background,tx,ty,fx,fy,[],Stacks)
%Stacks = a structure whose fields are consist of all the multi-resolution stacks used by the RV1 model. To create all
%         multi-resolution stacks, use the function: multi_resolution_stacks. retina_V1_model will create all necessary
%         multi-resolution stacks if 'Stacks' is not specified and the input is []. 
%         For example: X = retina_V1_model(target,background,tx,ty,fx,fy,Parameters,[])

%FIELDS OF 'Parameters':
%sl = standard deviation of a 2D Gaussian that specifies the region over which local luminance is calcualted. This 
%     remains constant as a function of retinal eccentricity. sl = 120 is the suggested default. Range: sl >= 1
%kc = scalar that specifies the size of the center region of a ganglion cell receptive field. Range: kc >= 1
%ks = scalar that specifies the size of the surround region of a ganglion cell receptive field. Range: ks > kc
%wc = weight on center in a Difference of Gaussians model of ganglion cell receptive fields. Range: [0,1]
%kb = weight on combined masking power of both narrowband and broadband masking components. Range: kb >= 0
%wb = weight on masking power of narrowband component (relative to broadband component). Range: [0,1]
%ro = pooling exponent. Range: ro > 0.
%p0 = masking power of uniform background (this parameter should change for every group of subjects). Range: p0 > 0
%dp_crit = d' criterion. We suggest setting dp_crit = 1.8009, which corresponds to 0.5+0.5*(1-exp(-1)) = 81.61%

%FIELDS OF 'Stacks':   
%T = target Optics and Gaussian stack
%B = background Optics and Gaussian stack
%E = target envelope stack
%DT = target Difference of Gaussians stack
%DB = background Difference of Gaussians stack
%FN = narrowband filter stack 

%DIRECTORY OF KEY FUNCTIONS
%To create a ganglion cell mosaic: ganglion_cell_mosaic
%To create any or all multi-resolution stacks: multi_resolution_stacks
%To modify the modulation transfer function for the human eye: MTF_optics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFAULT PARAMETERS

%Default pixels per degree
ppd = 120;

%Default ganglion cell spacing at fovea
s0 = spacing_fn(0,0);

%Default values for 'Parameters'
if isempty(Parameters)       
    Parameters.sl = 120; %Range: sl >= 1
    Parameters.kc = 1; %Range: kc >= 1
    Parameters.ks = 10.1; %Range: ks > kc
    Parameters.wc = 0.53; %Range: [0,1]
    Parameters.kb = 24; %Range: kb >= 0
    Parameters.wb = 0.925; %Range: [0,1]
    Parameters.ro = 2.4; %Range: p0 > 0    
    Parameters.p0 = 0.0014; %0.0014 for ModelFest, 0.00045 for our yes-no experiment    
    Parameters.dp_crit = 1.8009; %corresponds to 0.5+0.5*(1-exp(-1)) = 81.61 percent
end
P = Parameters; %rename

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GANGLION CELL MOSAIC

%To create a ganglion cell mosaic, use the function: ganglion_cell_mosaic

%Load previously created ganglion cell mosaic
load GC_mosaic; SAMP = cast(fliplr(GC_mosaic),'double');
SAMPx = center(size(SAMP,1)); SAMPy = center(size(SAMP,2)); %coordinates of mosaic center

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MULTI-RESOLUTION STACKS 

%To create all multi-resolution stacks, use the function: multi_resolution_stacks

target = cast(target,'double'); background = cast(background,'double'); %convert to doubles

if isempty(Stacks)
    %Create multi-resolution stacks 
    MRS = multi_resolution_stacks(target, background, P);
    T = MRS.target_stack; %target Optics and Gaussian stack 
    B = MRS.background_stack; %background Optics and Gaussian stack
    E = MRS.target_envelope_stack; %target envelope stack
    DT = MRS.DoG_target_stack; %target Difference of Gaussians stack
    DB = MRS.DoG_background_stack; %background Difference of Gaussians stack
    NF = MRS.narrowband_filter_stack; %narrowband filter stack 
else
    %Load and rename previously created multi-resolution stacks
    T=Stacks.target_stack; B=Stacks.background_stack; E=Stacks.target_envelope_stack; 
    DT=Stacks.DoG_target_stack; DB=Stacks.DoG_background_stack; NF=Stacks.narrowband_filter_stack; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOOP OVER TARGET AND FIXATION LOCATIONS

szB = size(B); cBx = center(szB(2)); cBy = center(szB(1)); %center of background   
szT = size(T); [BPx,BPy] = deal(2*szT(1)); %background patch is a square with twice the dimensions of the square target
cBPx = center(BPx); cBPy = center(BPy); %center coordinates of background patch

[C,C_db] = deal(zeros(1,length(fx)));
for qq = 1:length(fx)           
    
    %Use portion of stack corresponding to background patch   
    B_patch = cut_stack(B, cBy-round(ty(qq)*ppd), cBx+round(tx(qq)*ppd), BPy, BPx);
    DB_patch = cut_stack(DB, cBy-round(ty(qq)*ppd), cBx+round(tx(qq)*ppd), BPy, BPx);   
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FOVEATION AND GANGLION CELL SAMPLING
    
    %Target location relative to fixation
    ftx = fx(qq)-tx(qq); fty = fy(qq)-ty(qq);
    
    %Create resolution map at target location (right eye, nasal to the left)
    resmap = create_resmap(BPx, BPy, P.kc*s0, ftx, fty);
    resmap_cut = cut_stack(resmap, cBPx, cBPy, szT(1), szT(2));
    
    %Local luminance
    L = stack_interp(B_patch, P.sl);
    L = L(:,:,1); %size is fixed across eccentricity
    L_cut = cut_stack(L, cBPy, cBPx, szT(2), szT(1));
    
    %Foveated target    
    DT_fov = (resmap_interp(DT, resmap_cut))./L_cut;  
    
    %Foveated background    
    B_fov = (resmap_interp(B_patch, resmap))./L; 
    DB_fov = (resmap_interp(DB_patch, resmap))./L - (2.*P.wc-1); %subtract response to uniform background
    
    %Samplings function for target and background
    SAMP_cut_T = cut_stack(SAMP, SAMPy+round(fty*ppd), SAMPx+round(ftx*ppd), szT(2), szT(1));
    SAMP_cut_B = cut_stack(SAMP, SAMPy+round(fty*ppd), SAMPx+round(ftx*ppd), BPy, BPx);

    %Foveated envelope (resized to background)
    E_fov = resmap_interp(E, resmap_cut);
    E_fov = center_matrix(E_fov,zeros(BPx,BPy));    
    
    %Sampled envelope
    E_samp = E_fov.*SAMP_cut_B;    
    E_samp = E_samp./sum(E_samp(:));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DETECTION THRESHOLD
    
    %Pooled target responses
    T_pool = (sum(sum(abs(SAMP_cut_T.*DT_fov).^P.ro))).^(1./P.ro);
    
    %Broadband power
    BB = sum(sum(E_samp.*(DB_fov.^2)));
    
    %Narrowband power
    resmap_F = resmap(cBPx,cBPy)*ones(BPx, BPy); %resolution map for narrowband filter
    NF_ip = resmap_interp(NF, resmap_F); %Interpolated narrowband filter
    fft_BZ = fft2(ifftshift(B_fov)); 
    R_nb = abs(real(fftshift(ifft2(NF_ip.*fft_BZ)))); %Apply narrowband filter to background in Fourier space
    NB = sum(sum(E_samp.*(R_nb.^2))); %Narrowband power
    
    %Effective background noise    
    N_eff = sqrt(P.p0 + P.kb.*P.wb.*NB + P.kb.*(1-P.wb).*BB);
    
    %Predicted contrast threshold
    C(qq) = P.dp_crit.*N_eff./T_pool;
    C_db(qq) = 20*log10(C(qq)); %in decibels  
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT

out.target = target;
out.background = background;
out.tx = tx;
out.ty = ty;
out.fx = fx;
out.fy = fy;
out.Parameters = P;
out.threshold = C;
out.decibel_threshold = C_db;

