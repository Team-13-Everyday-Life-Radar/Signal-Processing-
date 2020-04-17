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

%2. Enable automatic trigger with frame time 1s
 oRS.oEPRadarBase.set_automatic_frame_trigger(5000000); 
 
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

minrange = 0.90; maxrange= 10;
minbeatfreq = ((chirpslope)*2*minrange)/c ;
maxbeatfreq = ((chirpslope)*2*maxrange)/c ;

threshhold = 25; % fft mag



% 
% figure (1)
%  plot(t,f)
%  grid on 
%  ylabel ('frequency(Hz)'); xlabel('time(s)'); 
%  title (' frequency of transmitted chirp');

Fs= 100000; %10kHz
[song, Fs] = audioread('Recording.mp3');
[song2, Fs] = audioread('Recording (2).m4a');
[song3, Fs] = audioread('Recording (3).m4a');
[song4, Fs] = audioread('Recording (4).m4a');
[song5, Fs] = audioread('Recording (5).m4a');
[song6, Fs] = audioread('Recording (6).m4a');
[song7, Fs] = audioread('Recording (7).m4a');
[song8, Fs] = audioread('Recording (8).m4a');
[song9, Fs] = audioread('Recording (9).m4a');
[song10, Fs] = audioread('Recording (10).m4a');
[song11, Fs] = audioread('Recording (11).m4a');
[song12, Fs] = audioread('Recording (12).m4a');
[song13, Fs] = audioread('Recording (13).m4a');
[song14, Fs] = audioread('Recording (14).m4a');
[song15, Fs] = audioread('Recording (15).m4a');
[song16, Fs] = audioread('Recording (16).m4a');
[song17, Fs] = audioread('Recording (17).m4a');
[song18, Fs] = audioread('Recording (18).m4a');


player = audioplayer(song,Fs);
player2 = audioplayer(song2, Fs);
player3 = audioplayer(song3, Fs);
player4 = audioplayer(song4, Fs);
player5 = audioplayer(song5, Fs);
player6 = audioplayer(song6, Fs);
player7 = audioplayer(song7, Fs);
player8 = audioplayer(song8, Fs);
player9 = audioplayer(song9, Fs);
player10 = audioplayer(song10, Fs);
player11 = audioplayer(song11, Fs);
player12 = audioplayer(song12, Fs);
player13 = audioplayer(song13, Fs);
player14 = audioplayer(song14, Fs);
player15 = audioplayer(song15, Fs);
player16 = audioplayer(song16, Fs);
player17 = audioplayer(song17, Fs);
player18 = audioplayer(song18, Fs);




while true
    % 3. Trigger radar chirp and get the raw data
    [mxRawData, sInfo] = oRS.oEPRadarBase.get_frame_data;
    ydata = mxRawData; % get raw data
  
    Amp(:,1,(1:cnum))= sqrt(real(ydata(:,1,(1:cnum))).^2 + imag(ydata(:,1,(1:cnum))).^2); 
    Amp2(:,1,(1:cnum))= sqrt(real(ydata(:,2,(1:cnum))).^2 + imag(ydata(:,2,(1:cnum))).^2); 

    %angle= atan(imag(ydata(:,1))./ real(ydata(:,1)));
    
    nfft = 2048;
    Af(:,1,(1:cnum)) = abs(fftshift(fft(Amp(:,1,(1:cnum)),nfft)));
    Af2(:,1,(1:cnum)) = abs(fftshift(fft(Amp2(:,1,(1:cnum)),nfft)));
    
    freq = (fs/nfft)*(-nfft/2:nfft/2 -1);
    fr = freq((nfft/2 +1):(end));
    
   y1 = mean(Af,3);
   y2 = mean(Af2,3);
   y= 0.5*(y1+y2);
  
%    figure(2)
%    plot(freq,y)
%     grid on 
%     ylabel ('amplitude'); xlabel('frequency (Hz)'); 
%     title (' spectrum of recevied signal');

   yhalf = y((nfft/2 +1):end);
   
%    figure(3)
%    plot(fr,yhalf)
%      grid on 
%     ylabel ('amplitude'); xlabel('frequency (Hz)'); 
%     title ('  one-side spectrum of recevied signal');

       stepfreq = fr(2) - fr(1); % affected by fs and nfft
       minindx = round((minbeatfreq - fr(1))/stepfreq);
       maxindx = round((maxbeatfreq - fr(1))/stepfreq); % rounded up or down ?? 
       
    frequ = fr((minindx+1):(maxindx+1));
    yfinal = yhalf((minindx+1):(maxindx+1));
%     
%     fixed y-axis needed
%     figure (4) 
%    plot(frequ,yfinal) 
%     grid on 
% ylabel ('amplitude'); xlabel('frequency (Hz)'); 
% title (' zoomed-in spectrum of the received signal');
%      axis = ([ frequ(1) frequ(end) 0 70]);
   

 [p , ind] = findpeaks(yfinal); %pass the averaged spectrum of all chirps
    
% p = p.*(p >= threshhold);
 
    for i= 1 : length(p)
        if p(i) >= threshhold 
            pks(i) = p(i);
            df(i)= frequ(ind(i)); %df(i) = (ind(i)-1).*(freq(2)-freq(1));
            target(i)= round((c*df(i))/(2*(chirpslope)),4,'significant');
           % target(i)= ceil(tgt(i)*4)/4;
        end
      
    end
    
    
     
   if (~exist('pks') || ~exist('df') || ~exist('target'))
           pks = 0;
          df = 0 ; 
          % target = 0;
        disp(' No Object Detected ')
        play(player)
        pause(3)
        stop(player)
   else
       target(target==0) = [];
       disp(round(target*4)/4)
        if round(target*4)/4 == 1
            play(player2)
            pause(3)
            stop(player2)
        elseif round(target*4)/4 == 1.25
            play(player3)
            pause(3)
            stop(player3)
        elseif round(target*4)/4 == 1.5
            play(player4)
            pause(3)
            stop(player4)
        elseif round(target*4)/4 == 1.75
            play(player5)
            pause(3)
            stop(player5)
        elseif round(target*4)/4 == 2
            play(player6)
            pause(3)
            stop(player6)
        elseif round(target*4)/4 == 2.25
            play(player7)
            pause(3)
            stop(player7)
        elseif round(target*4)/4 == 2.5
            play(player8)
            pause(3)
            stop(player8)
        elseif round(target*4)/4 == 2.75
            play(player9)
            pause(3)
            stop(player9)
        elseif round(target*4)/4 == 3
            play(player10)
            pause(3)
            stop(player10)
        elseif round(target*4)/4 == 3.25
            play(player11)
            pause(3)
            stop(player11)
        elseif round(target*4)/4 == 3.5
            play(player12)
            pause(3)
            stop(player12)
        elseif round(target*4)/4 == 3.75
            play(player13)
            pause(3)
            stop(player13)
        elseif round(target*4)/4 == 4
            play(player14)
            pause(3)
            stop(player14)
        elseif round(target*4)/4 ==4.25
            play(player15)
            pause(3)
            stop(player15)
        elseif round(target*4)/4 == 4.5
            play(player16)
            pause(3)
            stop(player16)
        elseif round(target*4)/4 == 4.75
            play(player17)
            pause(3)
            stop(player17)
        elseif round(target*4)/4 == 5
            play(player18)
            pause(3)
            stop(player18)
           
            
       end
   end
    
    figure (5)
    plot(frequ,yfinal,df,pks,'o')
    grid on
 ylabel ('amplitude'); xlabel('frequency (Hz)'); 
 title (' zoomed-in spectrum of the received signal with peaks circuled');
       axis = ([ frequ(1) frequ(end) 0 40]);

    rang = (c/(2*chirpslope)).*frequ;
    
    figure (6)
    plot(rang,yfinal)
    grid on
 ylabel ('amplitude'); xlabel('range(m)'); 
 title (' zoomed-in spectrum of the received signal');
    axis = ([ frequ(1) frequ(end) 0 40]);
%set(gca,'XTick',[rang(1):0.2:rang(end)]);%xicks=(0:0.5:10);


     clear pks; 
     clear target ;
    clear df ;

end

   
 
