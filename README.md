# mrQ_bb #
**mrQ** is a software package designed to calculate MR parameters (T1 and PD) using spoiled gradient echo scans (SPGR, FLASH). Using T1 and PD maps, mrQ performs the evaluation of macromolecular tissue volume (MTV) and the apparent volume of the interacting water protons (VIP) as well as the water-surface interaction rate (SIR). This is **version 2.1** of the mrQ software package.


The current repository has slight changes adapting th ecode to analyzing data of babies - which have smaller brains, and higher T1 values. 
It also includes the option of running the data with a low-resolution, unbiased T1 map (instead of using either the raw data or an aligned B1 map).

Please look at mrQ_baby.m for information on how to run a baby dataset 'traditionally'. 
Please look at mrQ_baby_seirT1map.m on how to run the data using a the SEIR-T1 map. s

To learn more about how to use mrQ please see the wiki:
https://github.com/mezera/mrQ/wiki

For more information, please contact: 

>Aviv Mezer: aviv.mezer(AT)elsc.huji.ac.il
>
>Shai Berman: shai.berman(AT)mail.huji.ac.il  

 
 
