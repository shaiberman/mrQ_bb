function [AnalysisInfo]=mrQ_SPM_reg_EPI2SPGR(AnalysisInfo,t1fileHM,outDir,B1file)
%[AnalysisInfo]=mrQ_ANTS_warp_EPI2SPGR(AnalysisInfo,t1fileHM,outDir,B1file)
%
% This function evaluates the ANTs affine linear registration from EPI to
% SPGR.  
%
% We apply the linear registration on the B1 map that was fitted in
% EPI space,  and bring it back to the SPGR space. The B1 will therefore be
% interpolated, but we are fine with that because we assume that B1 is a
% smooth parameter in the imaging space. The SPGR T1 is registered to the
% EPI warp version that was calculated before by
% mrQ_NLANTS_warp_SPGR2EPI_RB.m. This ensures that the registration is
% between two images with the same contrast. 
%
% The ANTs will register the EPI T1 to the SPGR T1 using unix functions.
% The code assumes ANTs and ITK code are part of the computer unix sell
% path (./bashrc). The registration happens in two steps: First, define the
% parameters; and second, apply it and warp the images. We warp both the T1
% and B1.
%
% The process is documented in the AnalysisInfo structure.
%
% INPUTS:
%
%      AnalysisInfo:  A structure that keeps records of the T1 fit files.
%
%          t1fileHM:  The path to SPGR T1 file. 
%                         (This is the uncorrected fast fit T1 map) 
%
%            outDir:  The path to where files are read from and written to
%
%            B1file:  The path to the B1 file in EPI space 
%                         (This file was saved by mrQ_smooth_LR_B1.m)
%
% OUTPUTS:
%
%      AnalysisInfo:  The updated structure that keeps records of the T1 
%                          fit files. The ANTs registration files' name and 
%                          date will be recorded in it.
%
% SEE ALSO: mrQ_T1M0_Fit.m and mrQ_smooth_LR_B1.m
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel
%   2015
%
%
%% I. make input structure

% % % s(1) =makeStructFromNifti(t1fileHM,-2);

s(1) = makeStructFromNifti(AnalysisInfo.T1_spgr_epi,-2);
s(2) = makeStructFromNifti(B1file,-2);

%% II. Align
 
interp = 1; %  1 (trilinear) or 7 (b-spline)
ref = readFileNifti(t1fileHM);
mmPerVox = ref.pixdim; % Set the resolution if it was not entered in

[s1,xform] = relaxFirstAlignAll(s,ref,mmPerVox,true,interp);

%% save output
AnalysisInfo.WARP_EPI_SPGR=fullfile(outDir,'Reg_EPI_SPGR');

% T1
AnalysisInfo.T1_epi_spgr=fullfile(outDir,'Warp_T1_EPI2SPGR.nii.gz');
dtiWriteNiftiWrapper(s1(1).imData, xform,AnalysisInfo.T1_epi_spgr);

% B1
AnalysisInfo.B1_epi_spgr=fullfile(outDir,'Warp_B1_EPI2SPGR.nii.gz');
dtiWriteNiftiWrapper(s1(2).imData, xform,AnalysisInfo.B1_epi_spgr);

 