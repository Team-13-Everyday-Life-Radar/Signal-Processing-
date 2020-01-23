clear all; close all;


c = physconst('LightSpeed');

waveform = phased.FMCWWaveform('SampleRate', 1e10,'SweepBandwidth',1e6,'SweepTime', 20e-5)%'SampleRate', 1e8,'SweepBandwidth',1e6)


T_sp = waveform.SweepTime; % sweep time (x1)
B=  waveform.SweepBandwidth ;  % bandwidth (x2)
R_max= 5; % maximum distance (x3)  
sweep_slope = B./T_sp; %0.00000000001*T_sp; %(x4)
fb_max = 2.*(2.*R_max).*B./(c.*T_sp); % maximum beat freq
fs = waveform.SampleRate; %10*fb_max; % sampling freq (x5)
%sig = step(waveform);

plot(waveform)

%%
obj_disChosen = 1.4 
obj_speed= 0;

object = phased.RadarTarget();
%sig = step(object,sig);

obj_motion = phased.Platform('InitialPosition',[obj_disChosen;0;0],...
'Velocity',[obj_speed;0;0]);

channel = phased.FreeSpace('PropagationSpeed',c,...
    'TwoWayPropagation',true);%%
%sig = step(channel,sig,[0;0;0],[obj_disChosen;0;0],[0;0;0],[obj_speed;0;0]);

%%

transmitter = phased.Transmitter();
%sig = step (transmitter, sig);

receiver = phased.ReceiverPreamp();
%sig = step (receiver, sig);


%%

 specanalyzer = dsp.SpectrumAnalyzer('SampleRate',fs,'PlotAsTwoSidedSpectrum',true,...
    'Title','Spectrum for received and dechirped signal',...
    'ShowLegend',true)


rng(1536);
%Nsweep = 64;
xr = complex(zeros(fs*T_sp,waveform.NumSweeps));

for m = 1:waveform.NumSweeps
    % Update targets positions
     [tgt_pos,tgt_vel] = obj_motion(T_sp); %[obj_disChosen;0;0][obj_speed;0;0]
    
    % Transmit FMCW waveform
    sig = waveform();
    txsig = transmitter(sig);
   
    
    % Propagate the signal and reflect off the target
    txsig = channel(txsig,[0;0;0],tgt_pos,[0;0;0],tgt_vel);
    txsig = object(txsig);

    % Dechirp the received radar return
    txsig = receiver(txsig);
    dechirpsig = dechirp(txsig,sig);

    % Visualize the spectrum
     specanalyzer([txsig dechirpsig]);

    xr(:,m) = dechirpsig;
end

%%


Dn = fix(fs/(299*fb_max));
for m = size(xr,2):-1:1
    xr_d(:,m) = decimate(xr(:,m),Dn,'FIR');
end

fs_d = fs/Dn;
fb_rng = rootmusic(pulsint(xr_d,'coherent'),1,fs_d);

rng_est = beat2range(fb_rng,sweep_slope,c)
% obj_distance = round(rng_est)

% %% Velocity 

%  peak_loc = val2ind(rng_est,c/(fs_d*2));
%  fd = -rootmusic(xr_d(peak_loc,:),1,fs_d);
%  v_est = dop2speed(fd,lambda)/2





