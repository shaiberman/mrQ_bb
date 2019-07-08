% This script organizes the input for running mrQ.
% In this script, we assume the input include VFAs for SPGRs and several ITs for the SEIR scans.
%%
clear,clc

%% define data dirs
mainDir = '../bb07';
indir = fullfile(mainDir,'input'); % here I have the spgr and seir images
outdir = fullfile(mainDir,'output'); % the mrQ will be saved here
addpath(genpath(fullfile(fileparts(mainDir),'code')))


%% define the SPGR hdr info:
inputData_spgr.rawDir =indir;
% A list of nifti names  (a unique string from the names is enough)
spnames = dir(fullfile(indir,'20670*.nii.gz'));% the number is changed according to the scan
for ii=1:length(spnames)
    inputData_spgr.name{ii}=spnames(ii).name;
    nii = readFileNifti(fullfile(indir,spnames(ii).name));
    
    st=regexp(nii.descrip,'fa=')+3;
    ed=regexp(nii.descrip,';ec')-1;
    inputData_spgr.flipAngle(ii)=str2double(nii.descrip(st:ed)); % The flip angle of each nifti in the list (degree)
end

inputData_spgr.TR = 14 * ones(1,length(spnames));  % the TR of each nifti in the list (msec)\
inputData_spgr.TE = 3 * ones(1,length(spnames)); % The TE of each nifti in the list (msec)
inputData_spgr.fieldStrength=3 * ones(1,length(spnames)); % The  field strength of each nifti in the list (Tesla)

%% define the SEIR hdr info:

inputData_seir.rawDir=indir;
% A list of nifti names  (a unique string from the names is enough)
senames = dir(fullfile(indir,'*epi*.nii.gz'));
for ii=1:length(senames)
    inputData_seir.name{ii}=senames(ii).name;
    
    st=regexp(senames(ii).name,'epi_')+4;
    ed=regexp(senames(ii).name,'.nii')-1;
    inputData_seir.IT(ii) =str2double( senames(ii).name(st:ed)); % The inversion time of each nifti in the list (msec)
end
inputData_seir.TR = 3000 * ones(1,length(senames));  % the TR of each nifti in the list (msec)
inputData_seir.TE = 40 * ones(1,length(senames));    % The TE of each nifti in the list (msec).

%% runnn

mrQ_run_bb(indir,outdir,inputData_spgr,inputData_seir,[],{'ants_bm',1,'refIm',nan,'testR1_BM',1,'AntsThresh',0.5})

% mrQ_run_bb calls some functions differently, such that the brainmask will
% be calculated differently and the assumptions on T1 values in WM will be
% higher. 

% the additional input parameter are the following:
  %  ants_bm=1, will make sure that when mrQ tries to register the SPGR T1
  %    map it will use the brainmasked map - this is necessary since the SEIR
  %    T1 map is also akull stripped. 
  %  refIm=nan, will make sure the code does not automaticall aligns the data
  %    to acpc space. 1st of all, it helps the registration to seir that they
  %    are in the similar spaces, 2nd,the acpc teplate is not for babies. 
  % AntsThresh=0.5, this is a threshold for the test of the registration
  %   between SEIR and SPGR. We usually change it according to the data set.
  %   it looks like in your data,fromthe several scans I analyzed, 0.5 is an
  %   ok threshold. this can be played with. 
  % testR1_BM=1, this uses the R1 (rather than M0, which is the default) to create the brainmask,
  %   because in this data it creates a saignificantly better brain mask.
