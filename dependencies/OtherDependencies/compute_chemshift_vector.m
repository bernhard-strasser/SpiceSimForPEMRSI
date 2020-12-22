function chemshift_vector = compute_chemshift_vector(LarmorFreqOrPar,dwelltime,vecSize,LowestPPMPoint)
% In units of:
% water_frequency: Hz
% dwelltime: s
% vecSize: -


%% Preparations

if(~exist('LowestPPMPoint','var'))
    LowestPPMPoint = [];
end
if(isstruct(LarmorFreqOrPar))
    LarmorFreq = LarmorFreqOrPar.LarmorFreq;
    if(~exist('dwelltime','var') || isempty(dwelltime))
        dwelltime = LarmorFreqOrPar.Dwelltimes(1)/1E9;
    end
    if(~exist('vecSize','var') || isempty(vecSize))
        vecSize = LarmorFreqOrPar.vecSize(1);
    end
else
    LarmorFreq = LarmorFreqOrPar;
end

%% 
%T_total = (vecSize - 1) * dwelltime;  % Total measurement time; The first point measuring at time 0, second point measuring at time dwelltime, ... So vecSize'th point measures at time (vecSize - 1) * dwelltime

%LarmorFreq = water_frequency/(1+4.65*10^-6);
% water_frequency = LarmorFreq * (1 + 4.65 * 10^-6);

%step_chemshift = 10^6/(LarmorFreq * T_total);	

% step_chemical_shift = chemshift_1 - chemshift_2 = (freq1 - TMS_freq)/TMS_freq*10^6 - (freq2 - TMS_freq)/TMS_freq*10^6 =
% = (freq1-freq2)/TMS_freq*10^6 = 1/T_total * 10^6/TMS_freq
% where chemshift_1 and 2 and freq1 and 2 are to chemshifts/freq of adjoining points


% so we know now the difference in chem shift between 2 points. But what chem. shift is eg the 512'th point?
% the most left point is the highes frequency, that is freq_max = water_frequency + bandwidth/2 (bandwidth symmetric about water_frequency)
% This corresponds to a chem shift of
% max_chemshift = (freq_max - TMS_freq)/TMS_freq * 10^6  = (water_freq - TMS_freq + bandwidth/2) / TMS_freq * 10^6
% min_chemshift = (freq_min - TMS_freq)/TMS_freq * 10^6  = (water_freq - TMS_freq - bandwidth/2) / TMS_freq * 10^6

bandwidth_frequency = 1/dwelltime;
step_frequency = bandwidth_frequency / vecSize;
freq_vector = ((ceil(vecSize/2)-1):-1:-floor(vecSize/2)) * step_frequency;

chemshift_vector = 10^6*freq_vector / LarmorFreq + 4.65;

% max_chemshift = (water_frequency - LarmorFreq + bandwidth_frequency/2)/(LarmorFreq)*10^6;
% min_chemshift = (water_frequency - LarmorFreq - bandwidth_frequency/2)/(LarmorFreq)*10^6;
% step_chemshift = (max_chemshift - min_chemshift)/(vecSize-1);
% 
% % The chemichal shift is now from min_chemshift to max_chemshift with steps step_chemshift
% 
% chemshift_vector = max_chemshift:-step_chemshift:min_chemshift; % High chemshift means small point because the FID got conjugated so that the spectrum is flipped left/right

if(~isempty(LowestPPMPoint))
   chemshift_vector = chemshift_vector - (min(chemshift_vector)-LowestPPMPoint); 
end

