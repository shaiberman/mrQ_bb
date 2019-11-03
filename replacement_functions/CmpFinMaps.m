%% sanity check:
% this script registers the final T1 map to the SEIR T1 map
% then, it looks at the histogram of the 2 map and creates a map of their difference
clear,clc
%% load mrQ
    
    spgrT1  = '/ems/elsc-labs/mezer-a/shai.berman/Documents/Code/testing_pipelines/mrQ_baby_Kalanit/bb08/output_RB/OutPutFiles_1/BrainMaps/T1_map_Wlin.nii.gz';
    seirT1  = '/ems/elsc-labs/mezer-a/shai.berman/Documents/Code/testing_pipelines/mrQ_baby_Kalanit/bb08/output_RB_someRatios/OutPutFiles_1/BrainMaps/T1_map_Wlin.nii.gz';
    % warpmat = mrQ.Ants_Info.WARP_SPGR_EPI; % the warp file that was reated by ANTS.
    seg     = '/ems/elsc-labs/mezer-a/shai.berman/Documents/Code/testing_pipelines/mrQ_baby_Kalanit/bb08/output_RB_someRatios/SPGR_1/Align_0.875_0.875_1/brainMask.nii.gz';
    
    ref = readFileNifti(seirT1);
    mmPerVox = ref.pixdim; % Set the resolution if it was not entered in
    
    %% REGISTER THE RELEVANT FILES TO seir
    % create a directory to ase your output
    outdir = '/ems/elsc-labs/mezer-a/shai.berman/Documents/Code/testing_pipelines/mrQ_baby_Kalanit/bb08/Warp_vsRB';
    if ~exist(outdir ,'dir'),mkdir(outdir),end
    interp = 1; %  1 (trilinear) or 7 (b-spline)
    % register the T1 map and segmentation file to the SEIR space andresample
    % to alow a fair comparoson
    
    % T1 map
    inf = spgrT1;
    t1out = fullfile(outdir,'T1_Wlin_map_Warped2SEIR.nii.gz');
    
    s(1) = makeStructFromNifti(spgrT1,-2);
    [s1,xform] = relaxFirstAlignAll(s,ref,mmPerVox,true,interp);
    dtiWriteNiftiWrapper(s1(1).imData, xform,  t1out);
    
    %% load the registered file and compare them
    T1h   = readFileNifti(t1out);   xform = T1h.qto_xyz; T1h   = T1h.data;
    T1ref = readFileNifti(seirT1);  T1ref = T1ref.data; % the seir is in ms, the spgr in sec..
    seg   = readFileNifti(seg);  seg   = seg.data;
    bm = logical(seg);
    
    %% first look at the histograms, in a plausible range of values
    
    msk = bm & T1h>1 & T1h<3 & T1ref>1 & T1ref<3;
       
    % find the ratio between the two peaks
    
    %% now look at the difference and and ratio of the maps
    
    t1diff = T1h - T1ref;
    t1diff(~bm)=nan;
    t1diff(isinf(t1diff))=nan;
    showMontage(t1diff)
    colorbar, colormap hot
    caxis([-0.2 0.2])
    title('Difference')
    dtiWriteNiftiWrapper(t1diff,xform, fullfile(outdir,'Difference.nii.gz'))
    
    t1r = T1h./T1ref;
    t1r(~bm)=nan;
    t1r(isinf(t1r))=nan;
    showMontage(t1r)
    colorbar, colormap hot
     caxis([0.9 1.1])
    title('Ratio')
    
