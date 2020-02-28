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
 oRS.oEPRadarBase.set_automatic_frame_trigger(3000000); %3 sec
% % 
%  oRS.oEPRadarBase.stop_automatic_frame_trigger; % stop it to change values
% % % 
% % % %  properties of radar
%    oRS.oEPRadarFMCW.lower_frequency_kHz = 24025000; % lower FMCW frequency 
%    oRS.oEPRadarFMCW.upper_frequency_kHz = 24225000; % upper FMCW frequency 
%    oRS.oEPRadarBase.num_chirps_per_frame = 1;  %[1,128]
%    oRS.oEPRadarBase.num_samples_per_chirp = 128; % [32, 64, 128, 256] 
%    oRS.oEPRadarADCXMC.samplerate_Hz = 640000; % 640000 affects the chirp time
%    oRS.oEPRadarFMCW.direction = 'Up Only';

c = 3e8;
smnum = double(oRS.oEPRadarBase.num_samples_per_chirp);
fs = double(oRS.oEPRadarADCXMC.samplerate_Hz);  
bw = double(oRS.oEPRadarFMCW.bandwidth_per_second);
ts =  double(oRS.oEPRadarBase.chirp_duration_ns); % chirp duration(ns) 
t =  0 : 1/fs : (1e-9*ts)-(1/fs);%linspace(0,1e-9.*ts,128); %time axis from 0 to 200 microsec;
f1= double(oRS.oEPRadarFMCW.lower_frequency_kHz)*1e3; % min freq 
f2= double(oRS.oEPRadarFMCW.upper_frequency_kHz)*1e3; % max freq
f= f1: (f2-f1)/smnum: f2-((f2-f1)/smnum); 
chirpslope = (f2-f1)/(t(end));

minrange = 0.95; maxrange= 5;
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
  
    A= sqrt(real(ydata(:,1)).^2 + imag(ydata(:,1)).^2);    
    angle= atan(imag(ydata(:,1))./ real(ydata(:,1)));

    nfft = 2048;
    
    freq = (fs/nfft)*(-nfft/2:nfft/2 -1);
    y = fftshift(fft(A,nfft));
    
    fr = freq((nfft/2 +1):(end));
    yhalf = y((nfft/2 +1):end);
   
       stepfreq = fr(2) - fr(1); % affected by fs and nfft
       minindx = round((minbeatfreq - fr(1))/stepfreq);
       maxindx = round((maxbeatfreq - fr(1))/stepfreq);
       
    frequ = fr((minindx+1):(maxindx+1));
    yfinal = yhalf((minindx+1):(maxindx+1));
    
    figure (2) 
    plot(freq,abs(y)) %plot(t, chirpsignal)
     grid on 
 ylabel ('amplitude'); xlabel('frequency (Hz)'); 
 title (' spectrum of recevied signal');
     
   figure (3) 
    plot(fr,abs(yhalf)) %plot(t, chirpsignal)
     grid on 
 ylabel ('amplitude'); xlabel('frequency (Hz)'); 
 title (' one-side spectrum of received signal');
 
    figure (4) 
   plot(frequ,abs(yfinal)) %plot(t, chirpsignal)
    grid on 
ylabel ('amplitude'); xlabel('frequency (Hz)'); 
title (' zoomed-in spectrum of the received signal');
     
    [p , ind] = findpeaks(abs(yfinal));
    
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
        disp(target)
   end
    
    figure (5)
    plot(frequ,abs(yfinal),df,pks,'o')
    grid on
 ylabel ('amplitude'); xlabel('frequency (Hz)'); 
 title (' zoomed-in spectrum of the received signal with peaks circuled');
    
    rang = (c/(2*chirpslope)).*frequ;
    
    figure (6)
    plot(rang,abs(yfinal))
    grid on
 ylabel ('amplitude'); xlabel('range(m)'); 
 title (' zoomed-in spectrum of the received signal');
 
    clear pks; 
    clear target ;
    clear df ;

end

   
 
