function [AnalysisInfo, Res]=mrQ_SPM_reg_SPGR2EPI(SET1file,T1mov,flipAngles,outDir,AlignFile,AligndSPGR,AnalysisInfo)


%mrQ_ANTS_warp_SPGR2EPI(AnalysisInfo,SET1file,flipAngles,outDir,AlignFile,AnalysisInfo)
% DOC!!!!
% NOTE: This function is currently misnamed. It is *not* a rigid body (RB)
% transformation.
%
%  This function evaluates the ANTs affine linear registration from SPGR to
%EPI. The ANTs will register the SPGR T1 to the SEIR EPI T1 using unix
%functions. The code assumes ANTs and ITK code are part of the computer
%unix sell path (./bashrc). The registration happens in two steps: First,
%define parameters; and second, apply it and warp the images. We warp both
%the T1 and the different flip angle degree SPGR raw images (all need to be
%aligned first).
%
% The process is documented in the AnalysisInfo structure.
%
% INPUTS:
%
%   SET1file     - The path from the SEIR T1 file to the SPGR will be
%                  registered to this target.
%
%
%   flipAngles   - The raw data flip angles. (The raw data are saved as NIfTI
%                  so the scan parameters are not saved in that file,
%                  but it is needed here).
%
%   outDir       - The path to which files are read from and written to.
%
%   AlignFile    - The name of the saved output file.
%
%   AligndSPGR   - The SPGR-aligned NIfTI files names.
%
%
% OUTPUTS:
%
%   AnalysisInfo - A structure that keeps records of the T1 fit files. The
%                  the ANTs registration files' name and date will be
%                  recorded in it.
%
%   Res          - A structure of registered files, ready for B1 fit. This
%                  structure is also saved as the AlignFile.
%
% see also: mrQfit_T1M0_ver2
%
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2016
%% I. make structure from niftii
s1(1) =makeStructFromNifti(T1mov,-2);
for ii = 1:numel(AligndSPGR)
    %  Load data from niftis - reshape and permute data (nifti) if needed
    
    [s1(ii+1)]= makeStructFromNifti([AligndSPGR{ii} '.nii.gz'],-2);
end
s = s1;
clear s1

%% II. Align

interp = 1; %  1 (trilinear) or 7 (b-spline)
ref = readFileNifti(SET1file);
mmPerVox = ref.pixdim; % Set the resolution if it was not entered in

[s1,xform] = relaxAlignAll(s,ref,mmPerVox,true,interp);

%% save output


% WARP_image= fullfile(outDir,'MoveIm');
for d=1:length(flipAngles) %loof over  flip Angles raw images
    % the name of the file we will align
    AnalysisInfo.Raw_spgr_epi{d}=[AligndSPGR{d} '_RB_SPGRT2EPI.nii.gz'];
    dtiWriteNiftiWrapper(s1(d+1).imData, xform,AnalysisInfo.Raw_spgr_epi{d});
end


%% V. Save and document

%make a structure to work with for B1 fit
t1seir=readFileNifti(SET1file);
Res{1}.im=t1seir.data;
Res{1}.name='t1SEIRepi';
Res{1}.xform=t1seir.qto_xyz;
clear t1seir
t1spgr=readFileNifti(AnalysisInfo.T1_spgr_epi);
Res{2}.im=t1spgr.data;
Res{2}.name='t1SPGR_in_epi_space';

for i=1:length(flipAngles)
    im=readFileNifti(AnalysisInfo.Raw_spgr_epi{i});
    Res{i+2}.im=im.data;
    Res{i+2}.name=['align_rawFA' num2str(flipAngles(i))] ;
    
end;


%save the structure
save(AlignFile,'Res');


%document and done
AnalysisInfo.spgr2epi_Align_date=date;
disp('DONE with linear registration SPGR to EPI ')

