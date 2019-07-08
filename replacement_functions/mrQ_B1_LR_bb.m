function mrQ=mrQ_B1_LR_bb(mrQ)
%  function mrQ=mrQ_B1_LR_tmp(mrQ)
%
% This function controls the B1 fit.
% The function has two parts: 
% 1. Fitting each voxel.
%     The fit is not independent to each voxel. A local regression (LR) 
%     approach is taken. When each voxel is fitted, close surrounding voxels 
%     are also used for the fit.
% 
% 2. Joining all the voxels to the B1 map.
% 	  After joining the voxels to the fits, we check for smoothness.  Then
%     we extrapolate to voxels that are missing in EPI space and then to
%     those in SPGR space.
%
% The INPUT is a mrQ structure, and the OUTPUT is the updated mrQ
% structure.
%
% Edited by A.M 2016
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel
% 2015
%

%% I. Get ready to fit B1 and call to fit
if ~isfield(mrQ,'B1fit_done');
    mrQ.B1fit_done=0;
end

B1_epi_OutPutDir=mrQ.Ants_Info.SEIR_SPGR_Curent_AlignDirs{2}; % the directory cotaining all images in the SPGR space 
AlignNum=mrQ.Ants_Info.SEIR_SPGR_Curent_AlignNums(1); % The index for the folder containing the best SEIR T1 fit for the registration. 

if ( mrQ.B1fit_done==0)
    
 %%   Register high-resolution EPI image to low-resolution aligned T1 image
   
    SPGR_niiFile_FA=eval(['mrQ.' mrQ.Ants_Info.SPGRFieldName '.SPGR_niiFile_FA']);
    AligndSPGR=eval(['mrQ.' mrQ.Ants_Info.SPGRFieldName '.AligndSPGR']);
    TRs=eval(['mrQ.' mrQ.Ants_Info.SPGRFieldName '.SPGR_niiFile_TR']);
 

%mrQ_NLANTS_warp_SPGR2EPI_RB(AnalysisInfo,SET1file,t1fileHM,flipAngles,outDir,AlignFile)

if ~isfield(mrQ,'SPGR_EPI_align_done');
    mrQ.SPGR_EPI_align_done=0;
end

if ( mrQ.SPGR_EPI_align_done==0)
    
 
    mrQ.spgr2epiAlignFile=fullfile(mrQ.Ants_Info.SEIR_SPGR_Curent_AlignDirs{2},'SEIRepiSPGRAlign_Struct.mat');
    %[mrQ.Ants_Info]=mrQ_NLANTS_warp_SPGR2EPI_RB(mrQ.SEIR_epi_T1file, mrQ.T1_LFit_HM, mrQ.SPGR_niiFile_FA, mrQ.spgr_initDir, mrQ.spgr2epiAlignFile, mrQ.AligndSPGR);
     [mrQ.Ants_Info]=mrQ_ANTS_warp_SPGR2EPI(mrQ.SEIRfits{AlignNum}.SEIR_epi_T1file,SPGR_niiFile_FA,B1_epi_OutPutDir,mrQ.spgr2epiAlignFile,AligndSPGR,mrQ.Ants_Info);
    mrQ.SPGR_EPI_align_done=1;
    
    save(mrQ.name,'mrQ');
    fprintf('\n Alignment of EPI to T1  - done!              \n');
else
    fprintf(['\n Using alignment of EPI to T1, calculated on '    mrQ.Ants_Info.spgr2epi_Align_date           '\n']);
    
end   
    

%% We build a mask for the voxel we'd like to fit in EPI space.

% [mrQ.maskepi_File] = mrQ_B1FitMask(mrQ.Ants_Info.dirAnts,mrQ.spgr2epiAlignFile,mrQ.SEIR_epi_fitFile,mrQ.spgr_initDir);
[mrQ.maskepi_File] = mrQ_B1FitMask([],mrQ.spgr2epiAlignFile,mrQ.SEIRfits{AlignNum}.SEIR_epi_fitFile,B1_epi_OutPutDir,mrQ.Ants_Info.WARP_SPGR_EPI);
    %define the fit parameters
    
    mrQ.B1.logname=mrQ_LB1SPGR_SEIR_Params(mrQ.name,mrQ.outDir,mrQ.sub,mrQ.spgr2epiAlignFile,mrQ.maskepi_File,TRs,SPGR_niiFile_FA);
    % save mrQ so the sungrid can load it fully characterized
    save(mrQ.name,'mrQ');
    % call to fit (with or without grid)    
    fprintf('\n fitting B1...              \n');

    mrQ_fitB1LR_Call(mrQ.B1.logname,mrQ.SunGrid);
    
    %check that the fit is done (for SGE)  before you move on (while loop)
    mrQ. B1fit_done=mrQ_Gridcheck(mrQ.B1.logname,mrQ.SunGrid,3);
    
    save(mrQ.name,'mrQ');
    
end
    
%% II. Build the grid B1 fits

    if isfield(mrQ,'B1Build_done');
    else
        mrQ.B1Build_done=0;
    end
    
    if (mrQ.B1Build_done==0 && mrQ.B1fit_done==1)
        fprintf('Build the B1 map from local fits    \n ');

        % join the estimations
        mrQ.B1=mrQ_build_LR_B1(mrQ.B1.logname,B1_epi_OutPutDir);
        
   % Smoothe the B1 fit, clear outliers, and extrapolate for voxels that 
   % have no fits.
        mrQ=mrQ_smooth_LR_B1_bb(mrQ,B1_epi_OutPutDir);
        
        mrQ=mrQ_build_epi2SPGR_B1_bb(mrQ,mrQ.Ants_Info.SEIR_SPGR_Curent_AlignDirs{2},mrQ.LinFit.T1_LFit);
        mrQ.B1Build_done=1;
        save(mrQ.name,'mrQ');
        
    end
            
%% correct B1 by a scalar
%  This part was found to be helpful in our data
%  We noticed that even after the B1  correction, the peaks of T1 in WM are
%  not aligned. In the ext section we atttempt to realign them.
%  It might require further improvement.

%     linearly fit T1 to align the peaks of the SPGR with that of the SEIR
if ~isfield(mrQ,'B1correct_done');
    mrQ.B1correct_done=0;
end
if  mrQ.B1correct_done==0;
    
    mrQ=mrQ_T1_LinFit(mrQ,mrQ.B1.unCorrected);
    
    % correct the B1
    mrQ=mrQ_correctB1_bb(mrQ);
    mrQ.B1correct_done=1;
    save(mrQ.name,'mrQ');
end
  
mrQ.B1_done=1; 
save(mrQ.name,'mrQ'); 
    %% non-linear fit of  B1
    
%     %%% FIT B1 by lsq fit compare T1 SEIR(Res{1}.im) to the multi flip
%     % angle Res{3:end}.im % USE sge make the fit faster
%     
%     intM0= double(mean(M0(brainMask)));
%     flipAngles = [s(:).flipAngle];
%     tr = [s(:).TR];
%     
%     if(~all(tr == tr(1))), error('TR''s do not match!'); end
%     tr = tr(1);
%     
%     % Load or Create the tissuemask from all non-zero points
%     tisuuemaskFile=fullfile(outDir,['maskepiF.nii.gz']);
%     
%     if exist(tisuuemaskFile,'file') && ~clobber
%         tisuuemask_=readFileNifti(tisuuemaskFile);
%         tisuuemask_=logical(tisuuemask_.data);
%         
%     else
%         
%         tisuuemask_ =zeros(size(Res{1}.im));
%         
%         % Binarize the mask
%         tisuuemask_(find(Res{1}.im>10 & Res{3}.im>0))=1;
%         
%         % Create a logical array from the tissue mask to index the
%         % non-zero locations
%         tisuuemask_=logical(tisuuemask_);
%         
%         % Save the tissue mask
%         dtiWriteNiftiWrapper(single(tisuuemask_), SE_Xform, tisuuemaskFile);
%         
%     end;
%     
%     if clobber && (exist([outDir '/tmpSG'],'dir'))
%         % in the case we start over and there are old  fits we will
%         % deleat them
%         eval(['! rm -r ' outDir '/tmpSG']);
%     end
%     
%     % USE sge make the B1 fit faster
%     
%     
%     % This is lsq fit that uses the grid but you can make it not use
%     % SGE: see help inside mrQ_fitB1_LSQ
%     [B1 resNorm dd] = mrQ_fitB1_LSQ(Res, tisuuemask_, tr,flipAngles, outDir, intM0, SE_Xform, SunGrid, 1, [sub 'B1fit'],mrQ.proclus);
%     dtiWriteNiftiWrapper(single(B1), SE_Xform, B1epifile);
%     dtiWriteNiftiWrapper(single(resNorm), SE_Xform, B1epiResidfile);
%     
% %end
% 
% %% Move B1 from epi space to SPGR smooth and upsample
% 
% 
% 
% % Smooth the SEIR B1 map by local and then fill the gap by global
% % regressions and then register back the smoothed B1 to the SPGR space using ANTS
% [B1,AnalysisInfo]=mrQ_smmothL_B1(B1epifile,AlignFile,outDir,B1file,xform,[],AnalysisInfo,t1fileHM,SET1Fitfile,B1epiResidfile);
