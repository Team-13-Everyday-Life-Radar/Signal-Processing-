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
% % % %  properties of radar
%    oRS.oEPRadarFMCW.lower_frequency_kHz = 24025000; % lower FMCW frequency 
%    oRS.oEPRadarFMCW.upper_frequency_kHz = 24225000; % upper FMCW frequency 
%    oRS.oEPRadarBase.num_chirps_per_frame = 1;  %[1,128]
%    oRS.oEPRadarBase.num_samples_per_chirp = 128; % [32, 64, 128, 256] 
%    oRS.oEPRadarADCXMC.samplerate_Hz = 640000; % 640000 affects the chirp time
%    oRS.oEPRadarFMCW.direction = 'Up Only';

c = 3e8;
fs = double(oRS.oEPRadarADCXMC.samplerate_Hz);  
bw = double(oRS.oEPRadarFMCW.bandwidth_per_second);
ts =  double(oRS.oEPRadarBase.chirp_duration_ns); % chirp duration(ns) 
t =  0 : 1/fs : (1e-9*ts)-(1/fs);%linspace(0,1e-9.*ts,128); %time axis from 0 to 200 microsec;
%t1 =  0 : 1/(4*f2) : (1e-9*ts)-(1/4*f2);%linspace(0,1e-9.*ts,128); %time axis from 0 to 200 microsec;
f1= double(oRS.oEPRadarFMCW.lower_frequency_kHz)*1e3; % min freq 
f2= double(oRS.oEPRadarFMCW.upper_frequency_kHz)*1e3; % max freq
f= f1: (f2-f1)/128: f2-((f2-f1)/128); 
chirpslope = (f2-f1)/(t(end));
%chirpsignal = cos(2*pi*(0.5*(bw*1e6)*(t1.^2)+f1*t1));

% 
figure (1)
 plot(t,f)
 grid on 
 ylabel ('frequency(Hz)'); xlabel('time(s)'); 
 title (' frequency of transmitted chirp');
%  
% figure (2) 
%  plot(t, chirpsignal)
%  grid on
%  ylabel ('amplitude'); xlabel('time(s)'); 
%  title ('tranmitted chirp') ;

while true
    % 3. Trigger radar chirp and get the raw data
    [mxRawData, sInfo] = oRS.oEPRadarBase.get_frame_data;
    ydata = mxRawData; % get raw data
   % target = oRS.oEPTargetDetection.get_targets; 
     %disp((ydata));

     
    A= sqrt(real(ydata(:,1)).^2 + imag(ydata(:,1)).^2);    
    angle= atan(imag(ydata(:,1))./ real(ydata(:,1)));

    nfft = 2048;
    freq = (fs/nfft)*(-nfft/2:nfft/2 -1);
    %fq = 0 : 320 : fs-1;
    y = fftshift(fft(A,nfft));
    
    figure (2) 
    plot(freq,abs(y)) %plot(t, chirpsignal)
     grid on
     
    figure (3)
    plot(freq(nfft/2:end-1),abs(y(nfft/2:end-1)))
    grid on 
    
    yfinal = (abs(y((nfft/2)+1:(3*nfft)/4)) > 20).*(abs(y((nfft/2)+1:(3*nfft)/4)));
    nfreq = freq((nfft/2)+1:(3*nfft)/4);
    figure (4)
    plot(nfreq,yfinal)
    grid on 
    
    [pks , ind] = findpeaks(yfinal);
    
    df = (ind(:)-1).*(nfreq(2)-nfreq(1));
    
    for i = 1:length(df)
        target(i)= (c*df(i))/(2*(chirpslope));
    end
    
    disp(target);
    
    %IFreceived = A'.*cos(2*pi*(0.5*(bw/1e-3)*(t.^2)+f1*t)+ angle'); % ??
%     IFreceived2 = A'.*cos(2*pi.50000000.*t+ angle');


% 
% figure(3)
% plot (t,IFreceived)
% grid on
% ylabel ('amplitude'); xlabel('time(s)'); 

% 
% figure(4)
% plot (t,IFreceived2)
% grid on
% ylabel ('amplitude'); xlabel('time(s)');


% nfft = 1024; 
% v= chirpsignal;
% vf = (1/nfft)*fftshift(fft(v, nfft));
% y = (1/nfft)*fftshift(fft(IFreceived, nfft));
% % y= y(1:nfft/2);
% % pxx = y.*conj(y)/(nfft*nfft);
% freq = (fs/nfft)*(-nfft/2:nfft/2 -1);
% % 
%  figure (4)
%  plot(freq,abs(vf))
%  grid on 
%  
%  figure (5)
%  plot(freq, abs(y))
%  grid on  
% 


end

   
 



  %  r = real(ydata(:,1));
%     j = imag(ydata(:,1));
%     I = r'.*cos(2.*pi.*f.*t);
%     Q = j'.*sin(2.*pi.*f.*t);
