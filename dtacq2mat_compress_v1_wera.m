%% dtacq2mat_radcelf_compress_v1_wera.mat
%% used matWERA libraries
function dtacq2mat_compress_v1_wera(filein,plt)
% function filein=dtacq2mat_compress_v1_wera(filein,plt)
% converts dta file to RAW file and compressed mat file; change parameters in beginning to suit site
% v1 - 2016.07.05 - LRB
%
%
% 11/15/2019 Douglas Cahl -- updated for WERA .RAW file
%
% output for WERA .RAW is timechirp(IQ,antN,chirpN,datalength)
%
% ex. 12 antennnas, 2048 chirps with a 0.433333s chirp time the
% length of the decimated 370Hz is 160 for each chirp
% -> timechirp(1:2,1:12,1:2048,1:160)
%
% plt = 1 optional plot of RAW and decimated data at 2 sample chirps
if nargin < 2 % default no plot
    plt = 0;
end


% 
RC.c = 299702520.2626;  % speed of radio waves in moist tropical air
RC.IQ = 2;              % I and Q channels
RC.SHIFT_FRAC = 0.08;   % fractional overlap 

% compression factor
RC.COMP_FAC = 16;       % compression factor set here, 16 for WERA with Fs = 5907.7
RC.OVER = 1;            % first round of decimation -- unused here (uncomment lines 116 and 124 to enable)

%% Read WERA RAW file
disp(['reading ' filein ' ...'])
[WERA,t,tc,I,Q,Fs]=read_WERA_raw(filein); % Read WERA RAW file

%% copy header params
RC.NANT = WERA.NAnt_SORT;
file_chirps = WERA.SAMPZ; % for WERA
RC.samps_per_chirp = WERA.MT;
RC.dec_samps_per_chirp = RC.samps_per_chirp/(RC.COMP_FAC*RC.OVER);
MTC = RC.dec_samps_per_chirp;

RC.NCHAN = RC.NANT;
RC.gain = 1;
RC.NCHIRP = WERA.SAMPZ;
RC.MT = WERA.MT;
RC.Fc = WERA.FREQ;
RC.chirp = WERA.RATE; % for WERA
RC.RAN_OFF = WERA.RAN_OFF;
RC.RAN_OFF_m = WERA.RAN_OFF_m;

% 
SHIFT_POS=floor((RC.SHIFT_FRAC.*MTC)/2)*RC.OVER*RC.COMP_FAC;  %scale it up for full chirp
MTL=RC.samps_per_chirp + SHIFT_POS*2; %the length of the chirp, with the added samples from shiftpos after
MTCL=ceil(MTL/(RC.OVER*RC.COMP_FAC));  % the decimated length of the chirp with the added samples from shiftpos
SHIFT_POS_dec = (MTCL-MTC)/2;



%% time of RAW data start
tt = [str2double(filein(end-18:end-15)) ... % yyyy
      str2double(filein(end-14:end-12)) ... % ddd
      str2double(filein(end-11:end-10)) ... % HH
      str2double(filein(end-9:end-8))];     % MM

mtime = datenum(tt(1),1,0)+ tt(2) + tt(3)/24 +tt(4)/24/60; % matlab time

% output file name
fileout = [filein(1:end-4) '.mat'];
%% decimation
disp('starting decimation ...')
%%% read, manipulate, and write data
timechirp=int16(zeros(RC.IQ,MTC,RC.NANT,file_chirps));

for ichirp=1:file_chirps
    dataI = nan(MTC,RC.NANT);  % decimated chirp length 
    dataQ = nan(MTC,RC.NANT);
    for ichan=1:RC.NANT
        %%% window, decimate, and unapply window to chirp
        datai0      = squeeze(I(ichan,ichirp,:));
        dataq0      = squeeze(Q(ichan,ichirp,:));
        ochirpN     = RC.samps_per_chirp + 2*SHIFT_POS;
        odec_chirpN = RC.dec_samps_per_chirp + 2*SHIFT_POS_dec;
        starti      = SHIFT_POS_dec+1;
        endi        = odec_chirpN - SHIFT_POS_dec;
        if ichirp == 1  % no pre for first chirp
            datai_pre = [];
            dataq_pre = [];
            ochirpN     = RC.samps_per_chirp + SHIFT_POS;
            odec_chirpN = RC.dec_samps_per_chirp + SHIFT_POS_dec;
            starti      = 1;
            endi        = odec_chirpN - SHIFT_POS_dec;
        else
            datai_pre = squeeze(I(ichan,ichirp-1,end-SHIFT_POS+1:end));
            dataq_pre = squeeze(Q(ichan,ichirp-1,end-SHIFT_POS+1:end));
        end
        
        if ichirp == file_chirps % no post for last chirp
            datai_post = [];
            dataq_post = [];
            ochirpN     = RC.samps_per_chirp + SHIFT_POS;
            odec_chirpN = RC.dec_samps_per_chirp + SHIFT_POS_dec;
            starti      = SHIFT_POS_dec+1;
            endi        = odec_chirpN;
        else
            datai_post = squeeze(I(ichan,ichirp+1,1:SHIFT_POS));
            dataq_post = squeeze(Q(ichan,ichirp+1,1:SHIFT_POS));
        end
        datai = cat(1,datai_pre,datai0,datai_post);
        dataq = cat(1,dataq_pre,dataq0,dataq_post);
        
        w = blackmanharris(ochirpN);
        w1 = blackmanharris(odec_chirpN);
        
        % I data
        wc = datai(:).*w(:);
%         dc = decimate(wc,RC.OVER,'fir');        % first round
        ddc = decimate(wc,RC.COMP_FAC,'fir');   % second round
        ddc = ddc(:)./w1(:);
        dataI(:,ichan) = ddc(starti:endi);
        
          
        % Q data
        wc = dataq(:).*w(:);
%         dc = decimate(wc,RC.OVER,'fir');        % first round
        ddc = decimate(wc,RC.COMP_FAC,'fir');   % second round
        ddc = ddc(:)./w1(:);
        dataQ(:,ichan) = ddc(starti:endi);
    end
    timechirp(1,:,:,ichirp)=dataI; % I data for this chirp
    timechirp(2,:,:,ichirp)=dataQ; % Q data
end


%% save file
disp(['saving ' fileout ' ...'])
save(fileout,'timechirp','RC','mtime','file_chirps','WERA') % ,'data','datac')
disp('done')

%% optional plot
if plt == 1
    ant_plt = 1;
    chirp1 = 399;
    chirp2 = 400;

    tc = 1:RC.samps_per_chirp;
    tc = tc/max(tc); % normalize
    
    t2 = 1:RC.dec_samps_per_chirp;
    t2 = t2*RC.OVER*RC.COMP_FAC;
    t2 = t2/max(t2); % normalize
    
    figure
    % I
    subplot(211)
    plot([tc+chirp1 tc+chirp2],[squeeze(I(ant_plt,chirp1,:)) ; squeeze(I(ant_plt,chirp2,:))])
    hold on
    plot([t2+chirp1 t2+chirp2],[squeeze(timechirp(1,:,ant_plt,chirp1))  squeeze(timechirp(1,:,ant_plt,chirp2))]')
    legend('I_{RAW}','I_{decimated}')

    % Q
    subplot(212)
    plot([tc+chirp1 tc+chirp2],[squeeze(Q(ant_plt,chirp1,:)) ; squeeze(Q(ant_plt,chirp2,:))])
    hold on
    plot([t2+chirp1 t2+chirp2],[squeeze(timechirp(2,:,ant_plt,chirp1))  squeeze(timechirp(2,:,ant_plt,chirp2))]')
    legend('Q_{RAW}','Q_{decimated}')
    xlabel('chirp #')
end

end
