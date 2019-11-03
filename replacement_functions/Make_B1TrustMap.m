function [BW2,out_mask]=Make_B1TrustMap(analysis,sub_path,outPath)


dir_BM=dir(fullfile(analysis,sub_path,'SPGR_1'));
BM_path=fullfile(analysis,sub_path,'SPGR_1',dir_BM(3).name,'/brainMask.nii.gz');

% the B1_LR file
infile=fullfile(analysis,sub_path,'SPGR_1',dir_BM(3).name,'/B1_LR.nii.gz');

% The EPI-T1 wraped to SPGR
wrap=fullfile(analysis,sub_path,'SPGR_1',dir_BM(3).name,'/WARP_EPI_SPGRWarp.nii.gz');

% text file for the alignment
txt=fullfile(analysis,sub_path,'SPGR_1',dir_BM(3).name,'/WARP_EPI_SPGRAffine.txt');

% refrence for the resulution of the output
ref=fullfile(analysis,sub_path,'OutPutFiles_1/BrainMaps/T1_map_Wlin.nii.gz');

% the B1_LR_2SPGR output
if notDefined('outPath')
    outPath=fullfile3(analysis,sub_path,'SPGR_1',dir_BM(3).name);
end
outf=fullfile(outPath,'B1_LR_2SPGR.nii.gz');

interpMethod = '--use-BSpline';

%% Warping
cmWarp=['xterm -e WarpImageMultiTransform  3 ' infile  ' ' outf ' -R ' ref ' ' wrap ' ' txt ' ' interpMethod];
[status, ~] = system(cmWarp);

if status ~= 0
    cmWarp=['WarpImageMultiTransform  3 ' infile  ' ' outf ' -R ' ref ' ' wrap ' ' txt ' ' interpMethod];
    % Run the command in unix and get back status and results:
    [~, ~] = system(cmWarp,'-echo');
end

%% Fill holes 
b1=readFileNifti(outf);
BW = imbinarize(b1.data);
%se = strel('sphere',10);
%closeBW = imclose(single(BW),se);
closeBW=BW;
BW2 = imfill(closeBW,8,'holes');
out_mask=fullfile(analysis,sub_path,'mrQ/SPGR_1/',dir_BM(3).name,'/B1_LR_2SPGR_mask.nii.gz');
dtiWriteNiftiWrapper(double(BW2),b1.qto_xyz,out_mask);

end