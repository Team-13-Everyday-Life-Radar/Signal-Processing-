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
 oRS.oEPRadarBase.set_automatic_frame_trigger(3000000);
%  oRS.oEPRadarBase.stop_automatic_frame_trigger; % stop it to change values
% % % 
% % % %  
%    oRS.oEPRadarFMCW.lower_frequency_kHz = 24025000; % lower FMCW frequency 
%    oRS.oEPRadarFMCW.upper_frequency_kHz = 24225000; % upper FMCW frequency 
%    oRS.oEPRadarBase.num_chirps_per_frame = 1;  %[1,128]
%    oRS.oEPRadarBase.num_samples_per_chirp = 128; % [32, 64, 128, 256] 
%    oRS.oEPRadarADCXMC.samplerate_Hz = 640000;
%    oRS.oEPRadarFMCW.direction = 'Up Only';


fs = double(oRS.oEPRadarADCXMC.samplerate_Hz);  
bw = double(oRS.oEPRadarFMCW.bandwidth_per_second);
ts =  double(oRS.oEPRadarBase.chirp_duration_ns); % chirp duration(ns) 
t =  0 : 1/fs : (1e-9*ts)-(1/fs);%linspace(0,1e-9.*ts,128); %time axis from 0 to 200 microsec;
f1= double(oRS.oEPRadarFMCW.lower_frequency_kHz)*1e3;
f2= double(oRS.oEPRadarFMCW.upper_frequency_kHz)*1e3;
f =f1: (f2-f1)/128: f2-((f2-f1)/128); %linspace(double(oRS.oEPRadarFMCW.lower_frequency_kHz)*1e3,double(oRS.oEPRadarFMCW.upper_frequency_kHz)*1e3,128);
chirpslope = (f2-f1)/(t(end));
chirpsignal = cos(2*pi*(0.5*(bw/1e-3)*(t.^2)+f1*t));


figure (1)
 plot(t,f)
 grid on 
 ylabel ('frequency(Hz)'); xlabel('time(s)'); 
 title (' frequency of transmitted chirp');
 
figure (2) 
 plot(t, chirpsignal)
 grid on
 ylabel ('amplitude'); xlabel('time(s)'); 
 title ('tranmitted chirp') ;

while true
    % 3. Trigger radar chirp and get the raw data
    [mxRawData, sInfo] = oRS.oEPRadarBase.get_frame_data;
    ydata = mxRawData; % get raw data
%    target = oRS.oEPTargetDetection.get_targets; 
     disp((ydata));

     % % % % Write ydata to text file (code ONLY WORKS for one receiver a
     % time)
     frames = 1 % how many chirps per frame is the data set to?
     fileID = fopen('ydata.txt','w');
     for i = 1:frames
         fprintf(fileID,'%f %f\n',[real(ydata_write(:,1,i)),imag(ydata_write(:,1,i))].');
     end
     fclose(fileID);
     % % % % 
     
    A= sqrt(real(ydata(:,1)).^2 + imag(ydata(:,1)).^2);    
    angle= atan(imag(ydata(:,1))./ real(ydata(:,1)));
%     r = real(ydata(:,1));
%     j = imag(ydata(:,1));
%     I = r'.*cos(2.*pi.*f.*t);
%     Q = j'.*sin(2.*pi.*f.*t);
    IFreceived = A'.*cos(2*pi*(0.5*(bw/1e-3)*(t.^2)+f1*t)+ angle');


figure(3)
plot (t,IFreceived)
grid on
ylabel ('amplitude'); xlabel('time(s)'); 



% nfft = 1024; 
% y = fft(IFreceived, nfft);
% y= y(1:nfft/2);
% pxx = y.*conj(y)/(nfft*nfft);
% freq = fs*(0:nfft/2 -1)/nfft;%(0:1:length(IFreceived)-1)*fs/128;
% 
%  figure (4)
%  plot(freq,abs(pxx)/128)
% % 



end

   
 

