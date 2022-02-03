%Clear previous data
clf
clear
clc
%.mat file from WFDB toolbox as suggested by Professor
%TO LOAD ANY MIT-BIT FILE RUN BELOW COMMAND IN COMMAND WINDOW% Also
%include mcode folder in path for running that command
%wfdb2mat('mitdb/100')  %This will downlaod full length signal for loading.
MITBD= [2273 1865 2187 2084 2229 2572 2027 2137 1774 2532]; % SOURCE %https://archive.physionet.org/physiobank/database/html/mitdbdir/records.htm#100
EData=[];
Err=[];
for i=1:10
EData=[EData; num2str((i+99)) 'm' '.mat'];
OsignalECG1= load(EData(i,:));
OsignalECG2 = OsignalECG1.val(1,1:end); 
OsignalECG=OsignalECG2';
Fs = 360;                           % Fs is the sampling frequency%%%MIT-BIH signal is sampled at 360
T = 1/Fs;                           % Sampling rate or frequency
N = length(OsignalECG);
ti = [0:N-1]/Fs;
sam=(1:length(OsignalECG));
[b,a] = butter(6, [45 55]/(Fs/2), 'stop');    %Powerline interference removal 
OsignalECG = filtfilt(b, a, OsignalECG);      %Taught in practice
[b,a] = butter(6, [95 105]/(Fs/2), 'stop');
OsignalECG = filtfilt(b, a, OsignalECG);

baseline = medfilt1(OsignalECG, floor(0.2 * Fs));  %Baseline correction
baseline = medfilt1(baseline, floor(0.6 * Fs));    %Idea from practice
OsignalECG = OsignalECG - baseline;
win=8;
l= gausswin(2*win+1)./win;
ECG_data1=zeros(size(OsignalECG));             

for k=1:length(OsignalECG)                  %Signal smoothening
for i=-win:win
    if i>-k && i<(length(OsignalECG)-k+1)
        ECG_data1(k)= ECG_data1(k)+OsignalECG(k+i)*l(i+win+1);
    end
end
end
N1=length(ECG_data1);
t1 = (0 : N1-1) ;
      
%Analyzing Heart beat

[rpks,Rpos]= findpeaks(ECG_data1(1:end)); %I am checking peaks for half of the datset i.e 1801
Limit = 2*(rms(ECG_data1(1:end)));        %In other words as my input file was                                            %for 10 secs, Now it will be only for 5 secs
R=rpks>Limit;   

Rw=zeros(size(ECG_data1(1:end)));
Rw(Rpos(R))=max(ECG_data1(1:end));
Beats= sum(R);
Err= [Err (MITBD(i)-Beats)];         %Error in our signals 
HRate=(Beats/3600)*30;  
end
plot(sam,OsignalECG);
title ('Plot of the original ECG signal sample')
xlabel ('Samples')
ylabel ('Amplitute')
grid on;
figure;
plot(ECG_data1(1:end));
hold on
stem(Rw,'r','^')
title('R WAVES for full length for one full signal')
grid on
figure;
errorbar(MITBD, Err)  %Plotting original vs Calculated beats in Errorplot as suggested by the professor
title('Error Bar for Original vs Processed signal')
grid on