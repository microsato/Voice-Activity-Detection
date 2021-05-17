% Read in Sound Data:
clear;
[d,s]= audioread('C:\Users\micro\.QtWebEngineProcess\Desktop\DSP_Project\male_no_noise1.wav');
t = 0:1:length(d)-1;
t=t/s;
for i =1:1:length(t)-1
    if t(i)==10
        cut =i;
    end
end 
%for babble noise signals
% d = d(1:cut);
% [d1, s1] = audioread('C:\DSP Project\babble.wav');
% d1 = interp(d1, 6);
% d1 = 0.1*d1;
% l=length(d1);
% d1(l+1)=0;
% n=d1+d';
% d=n;

%for white noisy signals
d = d(1:cut);
SNR=10; %set SNR
z=awgn(d,SNR,'measured');
figure();
d=z;
figure()
d = d(1:cut);
t = t(1:cut);
plot(t,d)
title('Sound Data versus Time')
xlabel('Time[seconds]');
ylabel('Volts[V]');

%Organise The Sound Data into Frames:
frame_length = 0.01;            %10ms per frame
n_frames = 10/frame_length;     %no. of frames in sound data

sample_per_frame = frame_length*s; 
last_samp = 1;
for k=1:1:n_frames
   new_d(:,k)= d(last_samp:last_samp+sample_per_frame) ;
   last_samp = last_samp + sample_per_frame;
end

%new_d(1281,:)=[];              %uncomment if repetive terms
new_d(:,1)=[];                  %where every column is 1 frame (10ms long)
%check if first column is all 0

%Get Energy Per Frame:
s_d = new_d.^2;                 %square each term 
for i = 1:1:n_frames-1
   col = s_d(:,i);
   e(i) = sum(col)./sample_per_frame;   %get Energy of each frame(column)
end

%Calculate the Threshold and Classify Sound Data:
T =0;
Emax =0;                    %Initiate maximumm Energy=0
Emin =0;                    %Initiate minimum Energy=0
count =0;                   %counts number of non-active frames
num = 30;                   %number of frames

for i= 1:1:n_frames-1       %loop through each frame
    E = e(i);               
    if i == 1               %The first frame:
        T(i) = E;           %Set Threshold, Emin,Emax = Energy of Frame
        Emax = E;           
        Emin = E;
    else                    %for all other frames
        if E>Emax           
            Emax = E;       
        end
        if E<Emin
            Emin = E;
        end
        k = 0.95;
        T(i) = (1-k)*Emax + k*Emin;%Calculate Threshold
        if E>T(i)                  %If E exceeds Thresh
            VAD(i) = 0.005;          %Frame is Voice
            count = 0;             %First voice frame
        else                       %if E less Thresh
            if count<num           %if number of frames since voice not detected less than num
                VAD(i)= 0.005;       %Frame is voice
                count = count + 1; %increase voice frame counter
            else
                VAD(i)=0;          %Frame is noise
            end
        end
        Emin = Emin*1.001;         %Update Emin
    end 
end
%Plot the Energy and Voice Detection:
figure()
plot(e,'k')
hold on
plot(VAD,'r')
title('Energy and Speech Detection of Sound Data');
xlabel("Frames[n]");
legend("Energy of Signal","Speech Detection","Location" ,"southeast");

%Plot the Sound Data and Voice Detection:
VAD_time(1)=0;
for i = 1:1:length(VAD)-1
    VAD_time(i+1)= 0.01*i;
end
figure()
plot(t,d,'b')
hold on
plot(VAD_time,4.3*8.*VAD, 'r')
title('Speech Data Undergoing Speech-Detection');
xlabel("Time[s]");

legend("Sound Data","Speech Detection","Location" ,"southeast");
