load('ECGData.mat');
total=ECGData.Data;

Fs=128;
sample=1000;%
t=0:(1/Fs):(sample/Fs)-1/Fs;
subplot(3,1,1);
plot(t,arr(20,1:sample));
title("Arrytmias")
ylabel("mV")
subplot(3,1,2);
plot(t,chf(1,1:sample));
title("Congestive cardiac failure")
ylabel("mV")
subplot(3,1,3);
plot(t,nsr(1,1:sample))
title("Normal sinus rhythm")
xlabel("Time(secs)")
ylabel("mV")


%butterworth filter 
fc=1;% cut off frequency
fn=Fs/2; %nyquivst frequency = sample frequency/2;
order = 5; %5th order filter, high pass
[b, a]=butter(order,(fc/fn),'high');
%fvtool(b,a);

%notch filter
wo = 60/(Fs/2);  
bw = wo/35;
[b_n,a_n] = iirnotch(wo,bw);
%fvtool(b_n,a_n);
for i=1:size(ECGData.Labels(1))
    but_sig= filter(b,a,total(i,:));
    not_sig= filter(b_n,a_n,but_sig);
    total(i,:)=not_sig;
end
arr=ECGData.Data(1:96,:);
chf=ECGData.Data(97:126,:);
nsr=ECGData.Data(127:162,:);
labels1=ones(length(arr(:,1)),1);
labels2=ones(length(chf(:,1)),1)*2;
labels3=ones(length(nsr(:,1)),1)*3;
labels=[labels1; labels2 ;labels3];

%filter signal
% y = filter(b,a,x);
% ydft = fft(y);
% plot(abs(ydft))

%peak find
%figure
m=size(total);
features=zeros(m(1),16);
for k=1:m(1)
if k~=103 && k~=132 && k~=135 
disp(k);
y1=total(k,1:sample);
tm=0:(1/Fs):(length(y1)/Fs)-1/Fs;
%plot(tm,y1)
y1=y1-mean(y1);
y1=y1/max(abs(y1));
y = abs(y1).^2;
[qrspeaks,locs] = findpeaks(y,tm,'MinPeakHeight',0.5,...
    'MinPeakDistance',0.50);
mean_peaks=0.5*mean(qrspeaks);
for i=1:length(qrspeaks)
    if qrspeaks(i)<mean_peaks
        qrspeaks(i)=[];
        locs(i)=[];
    end
end

% figure
% plot(tm,y)
% hold on
% plot(locs,qrspeaks,'ro')
% xlabel('Seconds')
% title('R Peaks Localized by Wavelet Transform with Automatic Annotations')
diff=[];
for i=1:(length(locs)-1)
    diff(i)=locs(i+1)-locs(i);
end

if length(diff)>0
rr_interval_over=mean(diff);features(k,1)=rr_interval_over;%first feature(rr)
rr_var=std(diff);features(k,2)=rr_var;%2 variability
end
if length(locs)>2
rr_interval=locs(3)-locs(2);
end
if length(locs)==2
    rr_interval=locs(2)-locs(1);
end
if length(locs)==1
    rr_interval=0;
end
if length(locs)>2 && locs(2)>0
seg=0;%uint8(rr_interval*0.4*Fs);
segment=y1(locs(2)*Fs-seg:locs(3)*Fs-seg);
t=0:(1/(Fs)):(length(segment)/Fs)-1/Fs;
if length(segment)>20

%plot(t,segment);


s=std(segment);features(k,3)=s;%3
sk=skewness(segment);features(k,4)=sk;%4
kr=kurtosis(segment);features(k,5)=kr;%5
%r=rms(segment);%6
segmentLength=length(segment);
order=uint8(length(segment)/2);
[Pxx,F] = pmcov(segment,order,segmentLength,Fs);
%pburg(segment,order,Fs)
[peak_cov,loc_cov]=findpeaks(Pxx,'MinPeakProminence',0,'NPeaks',1,'SortStr','descend');%7,8
if (length(peak_cov))>0
  features(k,6)=peak_cov;
  features(k,7)=loc_cov;
end


n=2^nextpow2(length(segment));
y= hilbert(segment);
env = abs(y);
fft2 = abs(fft(env,n));
[peaks2,locs2] = findpeaks(fft2(1:n/2),'MinPeakProminence',0, 'NPeaks', 3,'SortStr', 'descend');
peak_1 = peaks2(1);features(k,8)=peak_1;%9
peak_2 = peaks2(2);features(k,9)=peak_2;%10
%peak_3 = peaks2(3);%11
loc1 = locs2(1);features(k,10)=loc1;%12
loc2 = locs2(2);features(k,11)=loc2;%13
%loc3 = locs2(3);%14


[imf,residual] = emd(segment);
fft_imf1 = abs(fft(abs(hilbert(imf(:,1))),[],n));
fft_imf2 = abs(fft(abs(hilbert(imf(:,2))),[],n));
fft_imf3 = abs(fft(abs(hilbert(imf(:,3))),[],n));
[peaks_emd1,locs_emd1] = findpeaks(fft_imf1,'MinPeakProminence',0, 'NPeaks', 1,'SortStr', 'descend');%15,16
features(k,12)=peaks_emd1;features(k,13)=locs_emd1;
%[peaks_emd2,locs_emd2] = findpeaks(fft_imf2,'MinPeakProminence',0, 'NPeaks', 1,'SortStr', 'descend');%17,18
[peaks_emd3,locs_emd3] = findpeaks(fft_imf3,'MinPeakProminence',0, 'NPeaks', 1,'SortStr', 'descend');%19,20
features(k,13)=peaks_emd3;features(k,14)=locs_emd3;
rms_noise=rms(residual);%21
std_noise=std(residual);
features(k,15)=rms_noise;
features(k,16)=std_noise;
end
end
end
end

lda=fitcdiscr(features,labels);
ldaClass=resubPredict(lda);
ldaREsuberr=resubLoss(lda);
[ldaResub,grpOrd]=confusionmat(labels,ldaClass);
Acc=(ldaResub(1,1)+ldaResub(2,2)+ldaResub(3,3))/(162);
disp(Acc);
%pasge 293
MdlConf=fitcecoc(features,labels);
scvClass=resubPredict(MdlConf);
[svmResubCm,grOrder]=confusionmat(labels,scvClass);
acc=(svmResubCm(1,1)+svmResubCm(2,2)+svmResubCm(3,3))/162;
