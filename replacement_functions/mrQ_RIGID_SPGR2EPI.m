function [WarpFiles]= mrQ_RIGID_SPGR2EPI(LowResIm,HighResIm,outDir,morefiles2_HR) 
%
%mrQ_NLANTS_SPGR2EPI(LowResIm,HighResIm,mask,flipAngles,outDir,AlignFile,AligndSPGR)
% !!!!!DOC
%
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
%   LowResIm     - The path to target for example the from the SEIR T1 file to the SPGR will be
%                  registered to this target.
%
%   HighResIm     - The path to high resulotion image. for example the SPGR T1 file (this is uncorrected fast-fit
%                  T1 map). This file is used to register to the SEIR target
%                  because they have similar contrast (T1 maps).

%   outDir       - The path to which files are read from and written to.
%
%
%   morefiles2_HR   - The path to the hi-resultion files names to warp
%

%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2015

%
%% I. make structure from niftii
interp = 1; %  1 (trilinear) or 7 (b-spline)

s   = makeStructFromNifti(HighResIm,-2);
ref = readFileNifti(LowResIm);
% checks and balances

% Set the resolution if it was not entered in
mmPerVox = ref.pixdim;

[s1,xform] = relaxAlignAll(s,ref,mmPerVox,true,interp); 

%% save output
 

% WARP_image= fullfile(outDir,'MoveIm');
for d=1:length(morefiles2_HR) %loof over  flip Angles raw images
    
    file=dir(morefiles2_HR{d});
    savefileN=fullfile(outDir,['Reg' file.name]);

    dtiWriteNiftiWrapper(s1(d).imData, xform,savefileN);
    
    WarpFiles{d}=savefileN;
end    


