clc
disp('******************************************************************');
addpath('C:\IFX_P2G-HW-SW_V1.0.2\Firmware_Software\Communication Library\ComLib_Matlab_Interface\RadarSystemImplementation'); % add Matlab API
clear all; %#ok<CLSCR>
close all
resetRS; % close and delete ports

% 1. Create radar system object
szPort = findRSPort;  %scan all available ports
oRS = RadarSystem(szPort); % setup object and connect to board

disp('Connected RadarSystem:');
oRS %#ok<*NOPTS>

% 2. Enable automatic trigger with frame time 1s
 oRS.oEPRadarBase.set_automatic_frame_trigger(1000000); %1 sec
%  oRS.oEPRadarBase.stop_automatic_frame_trigger; % stop it to change values
% % 
% % %  properties of radar
%    oRS.oEPRadarFMCW.lower_frequency_kHz = 24025000; % lower FMCW frequency 
%    oRS.oEPRadarFMCW.upper_frequency_kHz = 24225000; % upper FMCW frequency 
%    oRS.oEPRadarBase.num_chirps_per_frame = 16;  %[1,128]
%    oRS.oEPRadarBase.num_samples_per_chirp = 128; % [32, 64, 128, 256] 
%    oRS.oEPRadarADCXMC.samplerate_Hz = 640000; % 640000 affects the chirp time
%    oRS.oEPRadarFMCW.direction = 'Up Only';

c = 3e8;
snum = double(oRS.oEPRadarBase.num_samples_per_chirp); % sample per chirp
cnum = double(oRS.oEPRadarBase.num_chirps_per_frame); % chirp per frame
fs = double(oRS.oEPRadarADCXMC.samplerate_Hz);  % samppling freq
bw = double(oRS.oEPRadarFMCW.bandwidth_per_second); % chcirpslope in microsec
ts =  double(oRS.oEPRadarBase.chirp_duration_ns); % chirp duration(ns) 
t =  0 : 1/fs : (1e-9*ts)-(1/fs); 
f1= double(oRS.oEPRadarFMCW.lower_frequency_kHz)*1e3; % min freq 
f2= double(oRS.oEPRadarFMCW.upper_frequency_kHz)*1e3; % max freq
f= f1: (f2-f1)/snum: f2-((f2-f1)/snum); 
chirpslope = (f2-f1)/(t(end));

minrange = 0.95; maxrange= 10;
minbeatfreq = ((chirpslope)*2*minrange)/c ;
maxbeatfreq = ((chirpslope)*2*maxrange)/c ;

threshhold = 25; % fft mag



% 
figure (1)
 plot(t,f)
 grid on 
 ylabel ('frequency(Hz)'); xlabel('time(s)'); 
 title (' frequency of transmitted chirp');


while true
    % 3. Trigger radar chirp and get the raw data
    [mxRawData, sInfo] = oRS.oEPRadarBase.get_frame_data;
    ydata = mxRawData; % get raw data
  
    Amp(:,1,(1:cnum))= sqrt(real(ydata(:,1,(1:cnum))).^2 + imag(ydata(:,1,(1:cnum))).^2); 
    %angle= atan(imag(ydata(:,1))./ real(ydata(:,1)));
    
    nfft = 2048;
    Af(:,1,(1:cnum)) = abs(fftshift(fft(Amp(:,1,(1:cnum)),nfft)));
    
    freq = (fs/nfft)*(-nfft/2:nfft/2 -1);
    fr = freq((nfft/2 +1):(end));
    
   y = mean(Af,3);
  
   figure(2)
   plot(freq,y)
    grid on 
    ylabel ('amplitude'); xlabel('frequency (Hz)'); 
    title (' spectrum of recevied signal');

   yhalf = y((nfft/2 +1):end);
   
   figure(3)
   plot(fr,yhalf)
     grid on 
    ylabel ('amplitude'); xlabel('frequency (Hz)'); 
    title ('  one-side spectrum of recevied signal');

       stepfreq = fr(2) - fr(1); % affected by fs and nfft
       minindx = round((minbeatfreq - fr(1))/stepfreq);
       maxindx = round((maxbeatfreq - fr(1))/stepfreq); % rounded up or down ?? 
       
    frequ = fr((minindx+1):(maxindx+1));
    yfinal = yhalf((minindx+1):(maxindx+1));
    
    % fixed y-axis needed
    figure (4) 
   plot(frequ,yfinal) 
    grid on 
ylabel ('amplitude'); xlabel('frequency (Hz)'); 
title (' zoomed-in spectrum of the received signal');
     axis = ([ frequ(1) frequ(end) 0 70]);
     
 [p , ind] = findpeaks(yfinal); %pass the averaged spectrum of all chirps
    
    for i= 1 : length(p)
        if (p(i) > threshhold)
            pks(i) = p(i);
            df(i)= frequ(ind(i)); %df(i) = (ind(i)-1).*(freq(2)-freq(1));
            target(i)= round((c*df(i))/(2*(chirpslope)),2,'significant');
        end
      
     end
    
   if (~exist('pks') || ~exist('df') || ~exist('target'))
           pks = 0;
           df = 0 ; 
           target = 0;
        disp(' No Object Detected ')
    else
        disp(ceil(target*4)/4)
   end
    
    figure (5)
    plot(frequ,yfinal,df,pks,'o')
    grid on
 ylabel ('amplitude'); xlabel('frequency (Hz)'); 
 title (' zoomed-in spectrum of the received signal with peaks circuled');
    
    rang = (c/(2*chirpslope)).*frequ;
    
    figure (6)
    plot(rang,yfinal)
    grid on
 ylabel ('amplitude'); xlabel('range(m)'); 
 title (' zoomed-in spectrum of the received signal');
 
    clear pks; 
    clear target ;
    clear df ;

end

   
 
