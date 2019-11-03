%% sanity check:
% this script registers the final T1 map to the SEIR T1 map
% then, it looks at the histogram of the 2 map and creates a map of their difference
clear,clc
%% load mrQ
figure
outNames = {'output'};%_RB_someRatios'};%{'output_RB'};%reg_smlMsk','output_RBreg_allRatios','output_RBreg_10smooth','output_RBreg_1smooth','output_RBreg','output_pickyB1'};
for ii = 1:length(outNames)
    mrQpath =  ['/ems/elsc-labs/mezer-a/shai.berman/Documents/Code/testing_pipelines/mrQ_baby_Kalanit/bb08/',outNames{ii},'/mrQ_params.mat'];
    
    load(mrQpath)
    
    %% define relevant files
    
    spgrT1  = mrQ.maps.T1path; % this is the final T map in the outputfiles/brainmaps directory
    if isfield(mrQ,'SEIR_epi_T1bmfile_fromR1')
        seirT1  = mrQ.SEIR_epi_T1bmfile_fromR1; % SEIR_epi_T1bmfile_fromR1;   % this should be the T1 map in the SEIR space that was used for the B1 fit. I take either the one with or without the brainmask
    else
        seirT1  = mrQ.SEIR_epi_T1file;
    end
    % warpmat = mrQ.Ants_Info.WARP_SPGR_EPI; % the warp file that was reated by ANTS.
    seg     = mrQ.T1w_tissue; % this is a segmentation file - if you have a better one it would be better to use it.
    fitmsk  = mrQ.maskepi_File; % this is the mask that defines which voxels were used for the B1 fit;. it is in SEIR space
    
    ref = readFileNifti(seirT1);
    mmPerVox = ref.pixdim; % Set the resolution if it was not entered in
    
    %% REGISTER THE RELEVANT FILES TO seir
    % create a directory to ase your output
    outdir = fullfile(mrQ.outDir,'sanChck');
    if ~exist(outdir ,'dir'),mkdir(outdir),end
    interp = 1; %  1 (trilinear) or 7 (b-spline)
    % register the T1 map and segmentation file to the SEIR space andresample
    % to alow a fair comparoson
    
    % T1 map
    inf = spgrT1;
    t1out = fullfile(outdir,'T1_Wlin_map_Warped2SEIR.nii.gz');
    segout = fullfile(outdir,'T1wtissueSeg_Warped2SEIR.nii.gz');
     if ~exist(t1out,'file') |  ~exist(segout,'file')
    % % %     cmWarp=['xterm -e WarpImageMultiTransform  3 ' inf  ' ' t1out ' -R '  seirT1 ' ' warpmat 'Warp.nii.gz ' warpmat 'Affine.txt'];
    % % %     [a,b] = system(cmWarp);
    % % %     if a~=0, disp(b),end
    
    s(1) = makeStructFromNifti(spgrT1,-2);
    s(2) = makeStructFromNifti(seg,-2);
    [s1,xform] = relaxFirstAlignAll(s,ref,mmPerVox,true,interp);
    dtiWriteNiftiWrapper(s1(1).imData, xform,  t1out);
    dtiWriteNiftiWrapper(round(s1(2).imData), xform,  segout);
    
    end
    % % % % seg file
    % % % inf = seg;
    % % % if ~exist(segout,'file')
    % % %     % % %     cmWarp=['xterm -e WarpImageMultiTransform  3 ' inf  ' ' segout ' -R '  seirT1 ' ' warpmat 'Warp.nii.gz ' warpmat 'Affine.txt --use-NN'];
    % % %     % % %     [a,b] = system(cmWarp);
    % % %     % % %     if a~=0, disp(b),end
    % % %
    % % %
    % % % end
    
    %% load the registered file and compare them
    T1h   = readFileNifti(t1out);   xform = T1h.qto_xyz; T1h   = T1h.data;
    T1ref = readFileNifti(seirT1);  T1ref = T1ref.data./1000; % the seir is in ms, the spgr in sec..
    seg   = readFileNifti(segout);  seg   = seg.data;
    fmsk  = readFileNifti(fitmsk);  fmsk  = fmsk.data;
    bm = logical(seg);
    
    %% first look at the histograms, in a plausible range of values
    
    msk = bm & T1h>1 & T1h<3 & T1ref>1 & T1ref<3;
    
    % % figure, hold on
    % % histogram(T1h(msk))
    % % histogram(T1ref(msk))
    % % legend('SPGR-T1','SEIR-T1')
    
    % % % % % % % % % % figure, hold on
    % % % % % % % % % %
    % % % % % % % % % % [spgr_Dens, spgr_Vals] = ksdensity(T1h(msk));
    % % % % % % % % % % mxT1_spgr = spgr_Vals(spgr_Dens==max(spgr_Dens));
    % % % % % % % % % % plot(spgr_Vals,spgr_Dens)
    % % % % % % % % % %
    % % % % % % % % % % [seir_Dens, seir_Vals] = ksdensity(T1ref(msk));
    % % % % % % % % % % mxT1_seir = ( seir_Vals(seir_Dens==max(seir_Dens)) ) ;
    % % % % % % % % % % plot(seir_Vals,seir_Dens)
    % % % % % % % % % %
    % % % % % % % % % % legend('SPGR-T1','SEIR-T1')
    
    % find the ratio between the two peaks
    
    %% now look at the difference and and ratio of the maps
    %
    % % % % % % t1diff = T1h - T1ref;
    % % % % % % t1diff(~bm)=nan;
    % % % % % % t1diff(isinf(t1diff))=nan;
    % % % % % % showMontage(t1diff)
    % % % % % % colorbar, colormap hot
    % % % % % % caxis([-0.5 0.5])
    % % % % % % title('Difference')
    % % % % % % dtiWriteNiftiWrapper(t1diff,xform, fullfile(outdir,'Difference.nii.gz'))
    %
    % % % % % % % t1r = T1h./T1ref;
    % % % % % % % t1r(~bm)=nan;
    % % % % % % % t1r(isinf(t1r))=nan;
    % % % % % % % showMontage(t1r)
    % % % % % % % colorbar, colormap hot
    % % % % % % %  caxis([0.8 1.2])
    % % % % % % % title('Ratio')
    % %
    %% now look at the difference and ratio of the maps only within the mask
    % (outside the mask we might expect a bit of bias due to extrapolation)
    %
    t1diff = T1h - T1ref;
    t1diff(~fmsk)=nan;
    t1diff(isinf(t1diff))=nan;
    % showMontage(t1diff)
    im = squeeze(t1diff(:,:,30));
   %  subplot(2,3,ii)
   imagesc(im)
    colorbar, colormap hot
    caxis([-0.5 0.5])
    title(['Difference - ' outNames{ii}])
    % % % % % % %
    % % % % % % % %
    % t1r = T1h./T1ref;
    % t1r(~fmsk)=nan;
    % t1r(isinf(t1r))=nan;
    % showMontage(t1r)
    % colorbar, colormap hot
    %  caxis([0.8 1.2])
    % title('Ratio')
    % % % % % % % %
    % showMontage(T1h),title('SPGR')
    % showMontage(T1ref),title('SEIR')
end
% suptitle('bm')
%% plot the ratio with the fit mask on top of it

% % threshold = 0.5;
% % GrayIm=mat2gray(t1diff,[-0.5 0.5]);
% % overIm=fmsk;
% % Incmap=hot(256);
% % clusterThresh=0;
% % ShowMOverlayImage(GrayIm,overIm,threshold,Incmap,clusterThresh )
% % title('ratio')
% %
% % threshold = 0.5;
% % GrayIm=mat2gray(t1r,[0.5 1.5]);
% % overIm=fmsk;
% % Incmap=hot(256);
% % clusterThresh=0;
% % ShowMOverlayImage(GrayIm,overIm,threshold,Incmap,clusterThresh )
% % title('difference')

%% plot the ratio with the trust mask on top of it
% make trust map
% analysis='/ems/elsc-labs/mezer-a/shai.berman/Documents/Code/testing_pipelines/mrQ_baby_Kalanit/bb07';
% sub_path='output_changedEPImask';
% [BW2,out_mask]=Make_B1TrustMap(analysis,sub_path,outdir);
% %
% % T1 = readFileNifti(spgrT1);
% %
% % threshold = 0.5;
% % GrayIm=mat2gray(T1.data,[1 2.5]);
% % overIm=BW2;
% % Incmap=hot(256);
% % clusterThresh=0;
% % ShowMOverlayImage(GrayIm,overIm,threshold,Incmap,clusterThresh )
% % title('T1 with B1 mask')
% %
t1diff = T1h - T1ref;
    t1diff(~bm)=nan;
    t1diff(isinf(t1diff))=nan;
    showMontage(t1diff)
    caxis([-0.5 0.5])
    
    colormap hot