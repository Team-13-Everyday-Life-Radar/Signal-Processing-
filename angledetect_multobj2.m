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
 
 oRS.oEPRadarBase.stop_automatic_frame_trigger; % stop it to change values
% 
% %  properties of radar
   oRS.oEPRadarFMCW.lower_frequency_kHz = 24025000; % lower FMCW frequency 
   oRS.oEPRadarFMCW.upper_frequency_kHz = 24225000; % upper FMCW frequency 
   oRS.oEPRadarBase.num_chirps_per_frame = 16;  %[1,128]
   oRS.oEPRadarBase.num_samples_per_chirp = 64; % [32, 64, 128, 256] 
   oEPTargetDetection.max_range_cm = 300;
   oRS.oEPRadarADCXMC.samplerate_Hz = 640000; % 640000 affects the chirp time
   oRS.oEPRadarFMCW.direction = 'Up Only';

   oRS.oEPRadarBase.set_automatic_frame_trigger(1000000); %1 sec
   
antenna_spacing = 6.22e-3; % in meters
c = 3e8;
snum = double(oRS.oEPRadarBase.num_samples_per_chirp); % sample per chirp
cnum = double(oRS.oEPRadarBase.num_chirps_per_frame); % chirp per frame
fs = double(oRS.oEPRadarADCXMC.samplerate_Hz);  % samppling freq
bw = double(oRS.oEPRadarFMCW.bandwidth_per_second); % chirpslope in microsec
ts =  double(oRS.oEPRadarBase.chirp_duration_ns); % chirp duration(ns) 
t =  0 : 1/fs : (1e-9*ts)-(1/fs); 
f1= double(oRS.oEPRadarFMCW.lower_frequency_kHz)*1e3; % min freq 
f2= double(oRS.oEPRadarFMCW.upper_frequency_kHz)*1e3; % max freq
f= f1: (f2-f1)/snum: f2-((f2-f1)/snum); 
chirpslope = (f2-f1)/(t(end));

fC = (f2+f1)/2;
   
lambda = c / fC;

minrange = 0.95; maxrange= 3;
minbeatfreq = ((chirpslope)*2*minrange)/c ;
maxbeatfreq = ((chirpslope)*2*maxrange)/c ;

threshhold = 20; % fft mag



while true
    % 3. Trigger radar chirp and get the raw data
    [mxRawData, sInfo] = oRS.oEPRadarBase.get_frame_data;
    ydata = mxRawData; % get raw data
  
%     Amp1(:,1,(1:cnum))= sqrt(real(ydata(:,1,(1:cnum))).^2 + imag(ydata(:,1,(1:cnum))).^2); 
%     Amp2(:,1,(1:cnum))= sqrt(real(ydata(:,2,(1:cnum))).^2 + imag(ydata(:,2,(1:cnum))).^2); 
    %angle= atan(imag(ydata(:,1))./ real(ydata(:,1)));
    
    Amp1(:,1,(1:cnum))= ydata(:,1,(1:cnum)); 
    Amp2(:,1,(1:cnum))= ydata(:,2,(1:cnum)); 
    nfft = 2048;
    
    Af1(:,1,(1:cnum)) = abs(fftshift(fft(Amp1(:,1,(1:cnum)),nfft)));
    Af2(:,1,(1:cnum)) = abs(fftshift(fft(Amp2(:,1,(1:cnum)),nfft)));
    
    Pf1(:,1,(1:cnum)) = (fftshift(fft(Amp1(:,1,(1:cnum)),nfft)));
    Pf2(:,1,(1:cnum)) = (fftshift(fft(Amp2(:,1,(1:cnum)),nfft)));
      
    freq = (fs/nfft)*(-nfft/2:nfft/2 -1);
    fr = freq((nfft/2 +1):(end));
    
   y1 = mean(Af1,3);
   y2 = mean(Af2,3);
   
   x1 = mean(Pf1,3);
   x2 = mean(Pf2,3);

   y1_half = y1((nfft/2 +1):end);
   y2_half = y2((nfft/2 +1):end);
   
   x1_half = x1((nfft/2 +1):end);
   x2_half = x2((nfft/2 +1):end);
   
       stepfreq = fr(2) - fr(1); % affected by fs and nfft
       minindx = round((minbeatfreq - fr(1))/stepfreq);
       maxindx = round((maxbeatfreq - fr(1))/stepfreq); % rounded up or down ?? 
       
    frequ = fr((minindx+1):(maxindx+1));
    y1_final = y1_half((minindx+1):(maxindx+1));
    y2_final = y2_half((minindx+1):(maxindx+1));
    
    x1_final = x1_half((minindx+1):(maxindx+1));
    x2_final = x2_half((minindx+1):(maxindx+1));
     
    [p1 , ind1] = findpeaks(y1_final); %pass the averaged spectrum of all chirps
 
    [~,idx] = sort(p1,'descend');
%     p1(idx(1)) amplitude of peak at index location
%     ind1(idx(1)) location of peak at index location
 
%d_phi for largest peak
            ang_rx_1_1 = angle(x1_final(ind1(idx(1))));
            ang_rx_2_1 = angle(x2_final(ind1(idx(1))));
    
            d_phi_1 = ang_rx_1_1 - ang_rx_2_1;
            
            if (d_phi_1 <= 0)
                d_phi_1 = d_phi_1 + 2*pi;
            end
             d_phi_1 = d_phi_1 - pi;

%d_phi for 2nd largest peak    

            ang_rx_1_2 = angle(x1_final(ind1(idx(2))));
            ang_rx_2_2 = angle(x2_final(ind1(idx(2))));
    
            d_phi_2 = ang_rx_1_2 - ang_rx_2_2;
            
            if (d_phi_2 <= 0)
                d_phi_2 = d_phi_2 + 2*pi;
            end
             d_phi_2 = d_phi_2 - pi;

            target_angle_1 = (asin((d_phi_1 * lambda) / (antenna_spacing * (2*pi)))); % AOA in radians

            target_angle_deg_1 = ((target_angle_1) * 180 / pi) % AOA in degrees

            target_angle_2 = (asin((d_phi_2 * lambda) / (antenna_spacing * (2*pi)))); % AOA in radians

            target_angle_deg_2 = ((target_angle_2) * 180 / pi) % AOA in degrees
    

end