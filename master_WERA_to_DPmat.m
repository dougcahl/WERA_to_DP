%% example of processing WERA .RAW and .SORT data to compatible .mat format
% for the WERA compatible fork of lera_DPd
%
% make sure the matWERA directory is your path
% (https://github.com/gvoulgaris0/matWERA)
% add the wera_raw_data path if not using full paths below
%
%
% Douglas Cahl 2019
% dcahl@geol.sc.edu

%% WERA RAW file to compressed mat file compatibile for LERA/WERA_DP code
clear
plt = 1; % plot two chirps 
filein = '20193221553_gtn.RAW';

% decimate .RAW WERA file
dtacq2mat_compress_v1_wera(filein,plt);

%% WERA SORT file to mat file compatibile for LERA/WERA_DP code
clear
filein = '20193221609_gtn.SORT';

[WERA,t,R,I,Q]=read_WERA_sort(filein);

% copy some stuff to RC
mtime = datenum([WERA.DATE WERA.TIME]);
RC.c = 299702520.2626;

RC.NANT = WERA.NAnt_SORT;
RC.NCHAN = RC.NANT;
RC.gain = 1;
RC.NCHIRP = WERA.SAMPZ;
RC.Fc = WERA.FREQ;
RC.chirp = WERA.RATE; 
RC.RAN_OFF = WERA.RAN_OFF;
RC.RAN_OFF_m = WERA.RAN_OFF_m;

% save file
disp(['saving ' filein(1:end-5) '.mat'])
save([filein(1:end-5) '.mat'],'RC','WERA','t','R','I','Q','mtime')


