% This script organizes the input for running mrQ.
% In this script, we assume the input include VFAs for SPGRs and Hua's
% output for the SEIR scans. 
%%
clear,clc

%% define data dirs
mainDir = '/ems/elsc-labs/mezer-a/shai.berman/Documents/Code/testing_pipelines/mrQ_baby_Kalanit/bb08';
indir = fullfile(mainDir,'input'); % here I have the spgr images, and an SEIR diuretory containing the ouptut of Hua's script
outdir = fullfile(mainDir,'output_RB_someRatios_lrgMask'); % the mrQ will be saved here
addpath(genpath(fullfile(fileparts(mainDir),'code')))

%% define the SPGR hdr info:
inputData_spgr.rawDir =indir;
% A list of nifti names  (a unique string from the names is enough)
spnames = dir(fullfile(indir,'20695*.nii.gz')); % the number is changed according to the scan
for ii=1:length(spnames)
    inputData_spgr.name{ii}=spnames(ii).name;
    nii = readFileNifti(fullfile(indir,spnames(ii).name));
    
    st=regexp(nii.descrip,'fa=')+3;
    ed=regexp(nii.descrip,';ec')-1;
    inputData_spgr.flipAngle(ii)=str2double(nii.descrip(st:ed)); % The flip angle of each nifti in the list (degree)
end

inputData_spgr.TR = 14 * ones(1,length(spnames));  % the TR of each nifti in the list (msec)
inputData_spgr.TE = 3 * ones(1,length(spnames)); % The TE of each nifti in the list (msec)
inputData_spgr.fieldStrength=3 * ones(1,length(spnames)); % The  field strength of each nifti in the list (Tesla)

%% define the SEIR hdr info:
% mrQ expects the fit to be saved a certain way.
% Therfore, I will save the differnt fit parameters in a structure that
% will be saved in the output directpry to overcome some of the code
% expectations. 
% organize the data into a structure
t1 = readFileNifti(fullfile(indir,'SEIR','20695_3_1_t1fit_t1.nii.gz'));
a  = readFileNifti(fullfile(indir,'SEIR','20695_3_1_t1fit_a.nii.gz'));
b  = readFileNifti(fullfile(indir,'SEIR','20695_3_1_t1fit_b.nii.gz'));
res = readFileNifti(fullfile(indir,'SEIR','20695_3_1_t1fit_res.nii.gz'));
bm  = readFileNifti(fullfile(indir,'SEIR','20695_3_1_t1fit_brain_mask.nii.gz'));

ll_T1(:,:,:,1)=t1.data;
ll_T1(:,:,:,2)=a.data;
ll_T1(:,:,:,3)=b.data;
ll_T1(:,:,:,4)=res.data;

mask = bm.data;

SEIRoutDir = fullfile(outdir,'SEIR');
if ~exist(SEIRoutDir,'dir'),mkdir(SEIRoutDir),end

T1fitFile = fullfile(SEIRoutDir,'T1_SEIR_T1_BM.mat');
save(T1fitFile,'ll_T1','mask');

% Here I create a link to the brainmask and the t1 map so that they will be
% saved in the output dorectory as well.

% The file t1_bmxform is the T1 map I resaved using the xform from the
% brainmask. This can simply be done with the following lines:
 % t = readFileNifti(t1path);
 % bm = readFileNifti(bmpath);
 % dtiWriteNiftiWrapper(t.data,bm.qto_xyz, NewPathToT1file)

inMap = fullfile(indir,'SEIR','t1_bmxform.nii.gz'); % this path depends on the name of the input file (20695_3_1_t1fit_t1.nii.gz)
mapLink = fullfile(SEIRoutDir,'T1_SEIR_T1_BM.nii.gz'); % this path needs to remain as is, because the code then looks for this specific path
if ~exist(mapLink)
    [a,b]= system(['ln -s ' inMap ' ' mapLink ]);
    if a~=0, disp(b),end
end

inBM = fullfile(indir,'SEIR','20695_3_1_t1fit_brain_mask.nii.gz');% this path depends on the name of the input file
bmLink = fullfile(SEIRoutDir,'T1_SEIR_BM.nii.gz'); % this path needs to remain as is, because the code then looks for this specific path
if ~exist(bmLink)
    [a,b]= system(['ln -s ' inBM ' ' bmLink ]);
    if a~=0, disp(b),end
end

%% runnn

mrQ_run_bb(indir,outdir,inputData_spgr,[],[],{'ants_bm',1,'refIm',nan,'AntsThresh',0.5, 'T1map_seir',T1fitFile})
% mrQ_run_bb calls some functions differently, such that the brainmask will
% be calculated diofferently and the assumptions on T1 values in WM will be
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
  % T1map_seir=FITFILE, this lets the code know to skip the T1 fit and look
  %   for the file containing the fit of the SEIR data.
  