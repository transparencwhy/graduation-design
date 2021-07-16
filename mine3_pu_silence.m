%谱分析捕获，加了自适应避障，无图 可求误码率，可研究同步误差
clc
clear all
% load FHMsignal&FHMsignalW;
er = [];
cc = [];
for iii = 1:1000
%------------初始条件设置------------
Fs = 2.5*10^6;%采样频率
Rs = 2.5*10^3;%码元传输速率
N = 560;%50码元，1码元1000点采样，必须是7000的倍数
N_de = 100;%解跳码元数
Len = N*Fs/Rs;%总点数长度
x = randi([0 1],1,N);%产生随机数56个作为数据码元
t = 0:1/Fs:(Len-1)/Fs;
menxian = 29;
%------------基带信号离散化------------
signal = [];
for k = 1:N
    
    if x(k) == 0
        signaltemp = -ones(1,1000);
    else
        signaltemp = ones(1,1000);
    end
    signal = [signal signaltemp];
end
t1 = 0:1*Rs/Fs:(Len-1)*Rs/Fs;
%------------------------------------------------------------------------------------------------------------------------
% figure(1)
% % subplot(1,1,1)
% plot(t1,signal);
% axis([0 50 -2 2]);xlabel('码元');ylabel('幅度');title('源信息序列')
%------------2FSK调制------------
f0 = 5*10^3;
f1 = 10*10^3;%2FSK两载波频率f0=5KHz,f1=10KHz
% T0 = 500;
% T1 = 250;
% s0 = gensig('sin',T0,N*Fs/Rs-1,1);
% s1 = gensig('sin',T1,N*Fs/Rs-1,1);
s0 = cos(2*pi*f0.*t);
s1 = cos(2*pi*f1.*t);
% s0 = rot90(s0);
% s1 = rot90(s1);
y0 = s0.*sign(-signal + 1);
y1 = s1.*sign(signal + 1);
FSK = y0 + y1;
%------------------------------------------------------------------------------------------------------------------------
% figure(2)
% subplot(2,1,1);
% plot(t,FSK);title('2FSK调制信号时域波形');xlabel('时间(秒)');ylabel('幅度');
% axis([0 0.002 -2 2]);
% FSKW = fft(FSK,Len);
% t1 = 0:Fs/(Len-1):Fs;
% subplot(2,1,2);
% plot(t1,abs(FSKW))
% title('2FSK频谱图');xlabel('频率(Hz)');ylabel('幅度');
% axis([0 15000 -inf +inf]);
%------------m序列生成------------
% an = [1 0 0];%3位寄存器
% len = length(an);
% L = 2^len-1;%m序列长度为7
% m = zeros(1,L)
% m(1)= an(3);
% for i = 2:L
%     an0=mod((an(1)+an(3)),2)
%     an(2:len)=an(1:len-1);%右移1位
%     an(1)=an0;
%     m(i) = an(3);
% end

%------------频率点初始化------------
FHPlocal = [313 453 523 243 383 173 103]*10^3;%本地载波与发射端载波相差中频2Khz
FHP_tiaoxu = [315 455 525 245 385 175 105]*10^3;%Fs/4=625KHz  2FSK第一零点带宽5+2*2.5=10KHz 频谱平移后占25KHz 跳频间隔为70KHz 足够
FHP_shunxu = [105 175 245 315 385 455 525]*10^3;%Fs/4=625KHz  2FSK第一零点带宽5+2*2.5=10KHz 频谱平移后占25KHz 跳频间隔为70KHz 足够
% a1 = 200;b1 = 350;
% a2 = 400; b2 = 550;
% a3 = 600; b3 = 750;
% a4 = 800; b4 = 950;
% a5 = 1000; b5 = 1150;
% a6 = 1200; b6 = 1350;
% a7 = 1400; b7 = 1550;
a_seq = [200 400 600 800 1000 1200 1400];
b_seq = [350 550 750 950 1150 1350 1550];
% aa1 = 1488; bb1 = 1641;
% aa2 = 1600; bb2 = 1750;
% aa3 = 1712; bb3 = 1862;
% new_aa = [aa1 aa2 aa3];
% new_bb = [bb1 bb2 bb3];
aa_seq = [1572 1768 1964];
bb_seq = [1722 1918 2114];

% c1 = 500;d1 = 670;
% c2 = 900;d2 = 1070;
% c3 = 1300; d3 = 1470;
% c4 = 1700; d4 = 1870;
% c5 = 2100; d5 = 2270;
% c6 = 2500; d6 = 2670;
% c7 = 2900; d7 = 3070;
c_seq = [500 900 1300 1700 2100 2500 2900];
d_seq = [670 1070 1470 1870 2270 2670 3070];

% cc1 = 3076; dd1 = 3246;
% cc2 = 3300; dd2 = 3470;
% cc3 = 3524; dd3 = 3694;
% new_cc = [cc1 cc2 cc3];
% new_dd = [dd1 dd2 dd3];
cc_seq = [3244 3636 4028];
dd_seq = [3414 3806 4198];
%------------部分频带阻塞干扰------------
noise_filter = fir1(512,[0.08,0.16]);%100K到200KHz
%------------高斯白噪声产生------------
gaosibaizaosheng = wgn(1,Len,10);
%------------单点干扰产生------------
t_noise = 0:1/Fs:(Len-1)/Fs;
s1 = sin(2*pi*FHP_shunxu(1)*t_noise);
s2 = sin(2*pi*FHP_shunxu(2)*t_noise);
s3 = sin(2*pi*FHP_shunxu(3)*t_noise);
s4 = sin(2*pi*FHP_shunxu(4)*t_noise);
s5 = sin(2*pi*FHP_shunxu(5)*t_noise);
s6 = sin(2*pi*FHP_shunxu(6)*t_noise);
s7 = sin(2*pi*FHP_shunxu(7)*t_noise);
%------------噪音、干扰选择------------
% interference_or_noise = s1;
ganrao =  zeros(1,Len);
%------------信道监测------------
FFTW = 0:Fs/(10000-1):Fs;
ganrao_clip = ganrao(1,1:10000);
ganrao_clipW = fft(ganrao_clip,10000);
plot(FFTW,abs(ganrao_clipW));
ganrao_clipW = abs(ganrao_clipW);
volume_jiance = [max(ganrao_clipW(400:450)) max(ganrao_clipW(680:730)) max(ganrao_clipW(960:1010)) max(ganrao_clipW(1240:1290)) max(ganrao_clipW(1520:1570)) max(ganrao_clipW(1800:1850)) max(ganrao_clipW(2080:2130))];
volume_jiance = volume_jiance>1500;
volume_beiyong = [595 665 735 805 875 945 1015]*10^3;%七个备用频率
volume_beiyong_zhongpin = [593 663 733 803 873 943 1013]*10^3;
for i = 1:7%自适应规避干扰频点。
    if volume_jiance(i) == 1
        b = find(FHP_shunxu == (35+70*i)*10^3);
        FHP_shunxu(b) = volume_beiyong(1);
        b = find(FHP_tiaoxu == (35+70*i)*10^3);
        FHP_tiaoxu(b) = volume_beiyong(1);
        FHPlocal(b) = volume_beiyong_zhongpin(1);
        volume_beiyong(1:6) = volume_beiyong(2:7);%左移
        volume_beiyong_zhongpin(1:6) = volume_beiyong_zhongpin(2:7);%左移
        a_seq(i) = aa_seq(1);
        aa_seq(1:2) = aa_seq(2:3);%左移
        b_seq(i) = bb_seq(1);
        bb_seq(1:2) = bb_seq(2:3);%左移
        c_seq(i) = cc_seq(1);
        cc_seq(1:2) = cc_seq(2:3);%左移
        d_seq(i) = dd_seq(1);
        dd_seq(1:2) = dd_seq(2:3);%左移
    end
end
%------------跳频载波对应------------
% FHP_shunxu = [105 175 245 315 385 455 525]*10^3;%Fs/4=625KHz  2FSK第一零点带宽5+2*2.5=10KHz 频谱平移后占25KHz 跳频间隔为70KHz 足够
t2 = 0:1/Fs:(7000-1)/Fs;
FHP1 = cos(2*pi.*FHP_shunxu(1).*t2);
FHP2 = cos(2*pi.*FHP_shunxu(2).*t2);
FHP3 = cos(2*pi.*FHP_shunxu(3).*t2);
FHP4 = cos(2*pi.*FHP_shunxu(4).*t2);
FHP5 = cos(2*pi.*FHP_shunxu(5).*t2);
FHP6 = cos(2*pi.*FHP_shunxu(6).*t2);
FHP7 = cos(2*pi.*FHP_shunxu(7).*t2);

% FHP = [205 275 345 415 485 555 625]*10^3;%Fs/4=625KHz  2FSK第一零点带宽5+2*2.5=10KHz 频谱平移后占25KHz 跳频间隔为70KHz 足够
% t2 = 0:1/Fs:(7000-1)/Fs;
% FHP1 =  (1/sqrt(2))*(cos(2*pi.*FHP(1).*t2)+sin(2*pi.*FHP(1).*t2));
% FHP2 =  (1/sqrt(2))*(cos(2*pi.*FHP(2).*t2)+sin(2*pi.*FHP(2).*t2));
% FHP3 = (1/sqrt(2))*(cos(2*pi.*FHP(3).*t2)+sin(2*pi.*FHP(3).*t2));
% FHP4 = (1/sqrt(2))*(cos(2*pi.*FHP(4).*t2)+sin(2*pi.*FHP(4).*t2));
% FHP5 = (1/sqrt(2))*(cos(2*pi.*FHP(5).*t2)+sin(2*pi.*FHP(5).*t2));
% FHP6 = (1/sqrt(2))*(cos(2*pi.*FHP(6).*t2)+sin(2*pi.*FHP(6).*t2));
% FHP7 = (1/sqrt(2))*(cos(2*pi.*FHP(7).*t2)+sin(2*pi.*FHP(7).*t2));
%------------信号跳频调制------------
carrier = [];
MAXCLOCK = N*Fs/Rs;%根据时间t进行选取，为50000
CLOCK = 0;
step = 7000; %共56个码元，每码元0.0004s，7个码元一跳频，1/0.0004/7≈357跳每秒
mm = 1;
an = [1 0 0];
len = length(an);
FHMsignal = [];
pattern = [];
carrier_total = [];
while CLOCK < MAXCLOCK
    CLOCK = CLOCK + step;
    register = [an(1),an(2),an(3)];
    decimal_reg = register*[2^2;2^1;2^0];
    
    
    switch(decimal_reg)
         case(1)
        pattern = [pattern FHP_shunxu(1) FHP_shunxu(1) FHP_shunxu(1) FHP_shunxu(1) FHP_shunxu(1) FHP_shunxu(1) FHP_shunxu(1)];
        carrier = FHP1;
        carrier_total = [carrier_total FHP1];
         case(2)
        carrier = FHP2;   
        pattern = [pattern FHP_shunxu(2) FHP_shunxu(2) FHP_shunxu(2) FHP_shunxu(2) FHP_shunxu(2) FHP_shunxu(2) FHP_shunxu(2)];
        carrier_total = [carrier_total FHP2];
         case(3)
        carrier = FHP3;
        pattern = [pattern FHP_shunxu(3) FHP_shunxu(3) FHP_shunxu(3) FHP_shunxu(3) FHP_shunxu(3) FHP_shunxu(3) FHP_shunxu(3)];
        carrier_total = [carrier_total FHP3];
         case(4)
        carrier = FHP4;
        pattern = [pattern FHP_shunxu(4) FHP_shunxu(4) FHP_shunxu(4) FHP_shunxu(4) FHP_shunxu(4) FHP_shunxu(4) FHP_shunxu(4)];
        carrier_total = [carrier_total FHP4];
         case(5)
        carrier = FHP5;  
        pattern = [pattern FHP_shunxu(5) FHP_shunxu(5) FHP_shunxu(5) FHP_shunxu(5) FHP_shunxu(5) FHP_shunxu(5) FHP_shunxu(5)];
        carrier_total = [carrier_total FHP5];
         case(6)
        carrier = FHP6;
        pattern = [pattern FHP_shunxu(6) FHP_shunxu(6) FHP_shunxu(6) FHP_shunxu(6) FHP_shunxu(6) FHP_shunxu(6) FHP_shunxu(6)];
        carrier_total = [carrier_total FHP6];
         case(7)
        carrier = FHP7;
        pattern = [pattern FHP_shunxu(7) FHP_shunxu(7) FHP_shunxu(7) FHP_shunxu(7) FHP_shunxu(7) FHP_shunxu(7) FHP_shunxu(7)];
        carrier_total = [carrier_total FHP7];
    end
    FHMsignal = [FHMsignal FSK(1,CLOCK-step+1:CLOCK).*carrier;];
    an0=mod((an(1)+an(3)),2);
    an(2:len)=an(1:len-1);%右移1位
    an(1)=an0;
end
%------------添加干扰------------
FHMsignal = FHMsignal + ganrao;
%------------绘制跳频图案------------ 
%------------------------------------------------------------------------------------------------------------------------
% figure(3)
% plot(pattern/10^3,'s','markerfacecolor','r','markersize',18); %'s'方块，'markerfacecolor'标记填充,'r'填充红色，'markersize'标记
% %大小
% ylabel('kHz');
% xlabel('码元');
% title('跳频图案');
% axis([1 56 0 700]);
% grid on;
%------------加噪声、干扰------------
% FHMsignal = FHMsignal + interference_or_noise;
% FHMsignal = awgn(FHMsignal,5);
% FHMsignal = FHMsignal + cos(2*pi*205*10^3*t);
%------------跳频信号时域频域分析------------
% yes = 0;
% if FSK.*carrier_total == FHMsignal
%     yes = 1;
% else 
%     yes = 0;
% end
%------------------------------------------------------------------------------------------------------------------------
% figure(4)
% FHMsignalfftW = fft(FHMsignal,length(FHMsignal));
% t1 = 0:Fs/(Len-1):Fs;
% subplot(2,1,1);plot(t,FHMsignal);xlabel('时间(s)');ylabel('幅度');title('跳频调制信号时域波形');
% axis([0 0.0032 -inf +inf]);
% subplot(2,1,2);plot(t1/10^3,abs(FHMsignalfftW));xlabel('频率（KHz）');ylabel('幅度');title('跳频调制信号频域波形');

%------------短时傅里叶变换分析跳频信号------------
% [S W T] = mf_spectrogram(FHMsignal,7000,6000,2.5*10^6);
% my_pcolor( W , T , S );

% spectrogram(FHMsignal,1000,500,1024,2.5*10^6,'yaxis')%处理信号，单次采样点数，相邻两次采样重叠点数，fft点数，采样频率，变换坐标轴
%------------初试条件再设置-------------
% Fs = 2.5*10^6;%采样频率
% Rs = 2.5*10^3;%码元传输速率
% N = 560;%50码元，1码元1000点采样
% Len = N*Fs/Rs;%总点数长度
% x = randi([0 1],1,N);%产生随机数56个作为数据码元
% t = 0:1/Fs:(Len-1)/Fs;
% t2 = 0:Fs/(1000-1):Fs;

%------------下变频变中频准备------------
t1 = 0:1/Fs:(1000-1)/Fs;
% FHPlocal = [313 453 523 243 383 173 103]*10^3;%本地载波与发射端载波相差中频2Khz
% % local_carrier = cos(2*pi*FHPlocal(2).*t1);
% % adddmax = 1000;
% % ii = 0;
% % for i = 1 : 500
% local_carrierW = fft(local_carrier,length(local_carrier));
% plot(t2,abs(local_carrierW));
%------------跳频信号与本地载波信号混频成中频信号------------
% % demode_signal = FHMsignal(7001-500+i:7001-500+i+999).*local_carrier;
% %------------滤波前时域------------
% figure(1);
% subplot(2,1,1);plot(demode_signal);title('滤波前时域');
% %------------滤波前频域------------
% demode_signalW = fft(demode_signal,length(demode_signal));
% subplot(2,1,2);plot(t2,abs(demode_signalW));title('滤波前频域');
% axis([0 10^5 -inf +inf]);
% %------------前1000点跳频信号频谱------------
% figure(2);
% FHMsignalW = fft(FHMsignal(1:1000),1000);
% plot(t2,abs(FHMsignalW));title('前1000点跳频信号频谱');
%------------混频信号过滤波器------------
b = fir1(48,0.0136);%125KHz以上60dB衰减
% % bpf_signal = filter(b,1,demode_signal);
% %------------滤波后时域------------
% figure(3);
% subplot(2,1,1);plot(bpf_signal);title('滤波后时域');
% %------------滤波后频域------------
% bpf_signalW = fft(bpf_signal,length(bpf_signal));
% subplot(2,1,2);plot(t2,abs(bpf_signalW));title('滤波后频域');
%------------平方累加能量检测------------
% % addd = 0;
% % for d = 1:1000
% %     addd = addd + (bpf_signal(d))^2;
% % end
% % if addd < adddmax
% %     adddmax = addd;
% %     ii = i;
% % end
% % end
% % adddmax
% % ii
%------------捕获预设定------------
Send_Clock = randi([1,49000]);
Initial_Send_Clock = Send_Clock;
MAXCLOCK = N*Fs/Rs;
tongbu_ok = 0;
tongbuchaoshi = 0;
buhuo_ok = 0;
xujing = 0;
a(1) = 0;
a(2) =0;
a(3) = 0;
an = [a(1) a(2) a(3)];
% Fs = 2.5*10^6;%采样频率
% Rs = 2.5*10^3;%码元传输速率
% N = 560;%560码元，1码元1000点采样
% N_de = 101;%解调码元长度，此处变的话本地载波应重新设置长度
% Len = N*Fs/Rs;%总点数长度
% x = randi([0 1],1,N);%产生随机数560个作为数据码元
% x = [0 0 0 0 0 0 0 1 1 1 1 1 1 1];
t = 0:1/Fs:(Len-1)/Fs;
W0 = 0:Fs/(14000-1):Fs;
W1 = 0:Fs/(7000-1):Fs;
% start = 42001;
out_of_signal = 0;
while (buhuo_ok == 0||tongbu_ok == 0)%进入大循环
    if out_of_signal ==1
        break
    end
%------------捕获逻辑------------
while(buhuo_ok == 0)
    if Send_Clock + 5*7000> MAXCLOCK
%         fprintf('数据耗尽无法捕获\n')
        out_of_signal = 1;
        break
    end
%------------时间捕获------------
signal_all = FHMsignal(1,Send_Clock:Send_Clock + 13999);
signal_first = FHMsignal(1,Send_Clock:Send_Clock + 6999);
signal_second = FHMsignal(1,Send_Clock+3500:Send_Clock + 3500+6999);
signal_third = FHMsignal(1,Send_Clock+7000:Send_Clock + 7000+6999);
FFT0 = fft(signal_all,length(signal_all));
FFT1 = fft(signal_first,length(signal_first));
FFT2 = fft(signal_second,length(signal_second));
FFT3 = fft(signal_third,length(signal_third));

% figure(11)
% subplot(4,1,1)
% plot(abs(FFT0));
% subplot(4,1,2)
% plot(abs(FFT1));
% subplot(4,1,3)
% plot(W1,abs(FFT2));
% subplot(4,1,4)
% plot(W1,abs(FFT3));

volum_14000 = [max(abs(FFT0(c_seq(1):d_seq(1)))) max(abs(FFT0(c_seq(2):d_seq(2)))) max(abs(FFT0(c_seq(3):d_seq(3)))) max(abs(FFT0(c_seq(4):d_seq(4)))) max(abs(FFT0(c_seq(5):d_seq(5)))) max(abs(FFT0(c_seq(6):d_seq(6)))) max(abs(FFT0(c_seq(7):d_seq(7)))) ];%捕获时间同步，用顺序排列便于判高低频填寄存器
volum_7000_first = [max(abs(FFT1(a_seq(1):b_seq(1)))) max(abs(FFT1(a_seq(2):b_seq(2)))) max(abs(FFT1(a_seq(3):b_seq(3)))) max(abs(FFT1(a_seq(4):b_seq(4)))) max(abs(FFT1(a_seq(5):b_seq(5)))) max(abs(FFT1(a_seq(6):b_seq(6)))) max(abs(FFT1(a_seq(7):b_seq(7))))];
volum_7000_second = [max(abs(FFT2(a_seq(1):b_seq(1)))) max(abs(FFT2(a_seq(2):b_seq(2)))) max(abs(FFT2(a_seq(3):b_seq(3)))) max(abs(FFT2(a_seq(4):b_seq(4)))) max(abs(FFT2(a_seq(5):b_seq(5)))) max(abs(FFT2(a_seq(6):b_seq(6)))) max(abs(FFT2(a_seq(7):b_seq(7))))];
volum_7000_third = [max(abs(FFT3(a_seq(1):b_seq(1)))) max(abs(FFT3(a_seq(2):b_seq(2)))) max(abs(FFT3(a_seq(3):b_seq(3)))) max(abs(FFT3(a_seq(4):b_seq(4)))) max(abs(FFT3(a_seq(5):b_seq(5)))) max(abs(FFT3(a_seq(6):b_seq(6)))) max(abs(FFT3(a_seq(7):b_seq(7))))];
volum_7000 = [volum_7000_first;volum_7000_second;volum_7000_third];
[m,index] = max(volum_14000);
[m,index1] = max(volum_7000(:,index));
% judge1 = index
Send_Clock = Send_Clock +(index1-1)*3500;
%------------频率捕获------------

zoom_or_zip = fix(index/4);
an(2:3)=an(1:3-1);
an (1) = zoom_or_zip;%填1

Send_Clock = Send_Clock + 7000;
FFTFFT2 = fft(FHMsignal(1,Send_Clock:Send_Clock + 6999));
% figure(2)
% plot(W1,abs(FFTFFT2))
volum_7000_repeat = [max(abs(FFTFFT2(a_seq(1):b_seq(1)))) max(abs(FFTFFT2(a_seq(2):b_seq(2)))) max(abs(FFTFFT2(a_seq(3):b_seq(3)))) max(abs(FFTFFT2(a_seq(4):b_seq(4)))) max(abs(FFTFFT2(a_seq(5):b_seq(5)))) max(abs(FFTFFT2(a_seq(6):b_seq(6)))) max(abs(FFTFFT2(a_seq(7):b_seq(7))))];
[m,index22] = max(volum_7000_repeat);
zoom_or_zip = fix(index22/4);
an(2:3)=an(1:3-1);
an (1) = zoom_or_zip;%填2

Send_Clock = Send_Clock + 7000;
FFTFFT3 = fft(FHMsignal(1,Send_Clock:Send_Clock + 6999));
% figure(3)
% plot(W1,abs(FFTFFT3))
volum_7000_repeat = [max(abs(FFTFFT3(a_seq(1):b_seq(1)))) max(abs(FFTFFT3(a_seq(2):b_seq(2)))) max(abs(FFTFFT3(a_seq(3):b_seq(3)))) max(abs(FFTFFT3(a_seq(4):b_seq(4)))) max(abs(FFTFFT3(a_seq(5):b_seq(5)))) max(abs(FFTFFT3(a_seq(6):b_seq(6)))) max(abs(FFTFFT3(a_seq(7):b_seq(7))))];
[m,index33] = max(volum_7000_repeat);
zoom_or_zip = fix(index33/4);
an(2:3)=an(1:2);
an (1) = zoom_or_zip;%填3，至此已填满an寄存器

% figure(4)
% W10 = 0:Fs/(49000-1):Fs;
% FFTFFT10 = fft(FHMsignal(1,1:49000),49000);
% plot(W10,abs(FFTFFT10))

% index
% index22
% index33
%     
%     if i == 1
% start = start + (index1-1)*3500
% end

%------------闭环跟跳确认------------
MAXCLOCK_pu = 7000*3;
CLOCK_pu = 0;
step_pu = 7000;
while(CLOCK_pu < MAXCLOCK_pu)
    CLOCK_pu = CLOCK_pu + step_pu;
    
    Send_Clock = Send_Clock + 7000;
    an0=mod((an(1)+an(3)),2);
    an(2:3)=an(1:3-1);%闭环右移1位
    an(1)=an0;
    decimal_an = an*[2^2;2^1;2^0];
    FFTFFT_check = fft(FHMsignal(1,Send_Clock:Send_Clock + 6999));
    volum_7000_check = [max(abs(FFTFFT_check(a_seq(1):b_seq(1)))) max(abs(FFTFFT_check(a_seq(2):b_seq(2)))) max(abs(FFTFFT_check(a_seq(3):b_seq(3)))) max(abs(FFTFFT_check(a_seq(4):b_seq(4)))) max(abs(FFTFFT_check(a_seq(5):b_seq(5)))) max(abs(FFTFFT_check(a_seq(6):b_seq(6)))) max(abs(FFTFFT_check(a_seq(7):b_seq(7))))];
    [m,index_check] = max(volum_7000_check);
    if(decimal_an ~= index_check)
%         fprintf('发生虚警\n')
        xujing = 1;
        break
    end
    if(CLOCK_pu == MAXCLOCK_pu)
        buhuo_ok = 1;
        volum_7000_transmit = [max(abs(FFTFFT_check(a_seq(4):b_seq(4)))) max(abs(FFTFFT_check(a_seq(6):b_seq(6)))) max(abs(FFTFFT_check(a_seq(7):b_seq(7)))) max(abs(FFTFFT_check(a_seq(3):b_seq(3)))) max(abs(FFTFFT_check(a_seq(5):b_seq(5)))) max(abs(FFTFFT_check(a_seq(2):b_seq(2)))) max(abs(FFTFFT_check(a_seq(1):b_seq(1))))];%跳序定judge2
        [m,index_transmit] = max(volum_7000_transmit);
        judge2 = index_transmit;   
    end
end
if (buhuo_ok == 1)
    break
end
end





% judge2
tb = 0;
tongbu_start = Send_Clock;
tongbu_end = Send_Clock + 6999;

if Send_Clock + 14000 > MAXCLOCK
%     fprintf('数据耗尽无法同步\n')
    break
end
while(buhuo_ok == 1 && tb < 100 && tongbu_ok == 0)%尝试进入跟踪环节
    
    if tongbu_end > Len
        Send_Clock = tongbu_end + 1;
        break
    end
    k1 = 3500;
%     before_signal = FHMsignal(1,tongbu_start:tongbu_end - k1);
%     after_signal = FHMsignal(1,tongbu_start + k1:tongbu_end);
    t3 = 0 : 1/Fs: (3500-1)/Fs;
    local_carrier_tongbu = cos(2*pi*(FHPlocal(judge2)).*t3);
        tb = tb + 1;
demode_signal1_tongbu = FHMsignal(1,tongbu_start - 2:tongbu_end - k1 - 2).*local_carrier_tongbu; 
demode_signal2_tongbu = FHMsignal(1,tongbu_start - 1:tongbu_end - k1 - 1).*local_carrier_tongbu; 
demode_signal3_tongbu = FHMsignal(1,tongbu_start:tongbu_end - k1).*local_carrier_tongbu; 
demode_signal4_tongbu = FHMsignal(1,tongbu_start + 1:tongbu_end - k1 + 1).*local_carrier_tongbu; 
demode_signal5_tongbu = FHMsignal(1,tongbu_start + 2:tongbu_end - k1 + 2).*local_carrier_tongbu; 
demode_signal6_tongbu = FHMsignal(1,tongbu_start + 3:tongbu_end - k1 + 3).*local_carrier_tongbu; 

bpf_signal1_tongbu =filter(b,1,demode_signal1_tongbu);
bpf_signal2_tongbu =filter(b,1,demode_signal2_tongbu);
bpf_signal3_tongbu = filter(b,1,demode_signal3_tongbu);
bpf_signal4_tongbu =filter(b,1,demode_signal4_tongbu);
bpf_signal5_tongbu =filter(b,1,demode_signal5_tongbu);
bpf_signal6_tongbu =filter(b,1,demode_signal6_tongbu);
Energy_Detect1 = 0;
 for ii = 1:3500
            Energy_Detect1 = Energy_Detect1 + ((bpf_signal1_tongbu(ii))^2 + (bpf_signal2_tongbu(ii))^2+(bpf_signal3_tongbu(ii))^2+(bpf_signal4_tongbu(ii))^2+(bpf_signal5_tongbu(ii))^2+(bpf_signal6_tongbu(ii))^2)/6;
 end
% demode_signal0_tongbu = FHMsignal(1,tongbu_start + k1 - 3:tongbu_end - 3).*local_carrier_tongbu; 
demode_signal1_tongbu = FHMsignal(1,tongbu_start + k1 - 2:tongbu_end - 2).*local_carrier_tongbu; 
demode_signal2_tongbu = FHMsignal(1,tongbu_start + k1 - 1:tongbu_end - 1).*local_carrier_tongbu; 
demode_signal3_tongbu = FHMsignal(1,tongbu_start + k1:tongbu_end).*local_carrier_tongbu; 
demode_signal4_tongbu = FHMsignal(1,tongbu_start + k1 + 1:tongbu_end + 1).*local_carrier_tongbu; 
demode_signal5_tongbu = FHMsignal(1,tongbu_start + k1 + 2:tongbu_end + 2).*local_carrier_tongbu; 
demode_signal6_tongbu = FHMsignal(1,tongbu_start + k1 + 3:tongbu_end + 3).*local_carrier_tongbu; 
% bpf_signal0_tongbu =filter(b,1,demode_signal1_tongbu);
bpf_signal1_tongbu =filter(b,1,demode_signal1_tongbu);
bpf_signal2_tongbu =filter(b,1,demode_signal2_tongbu);
bpf_signal3_tongbu = filter(b,1,demode_signal3_tongbu);
bpf_signal4_tongbu =filter(b,1,demode_signal4_tongbu);
bpf_signal5_tongbu =filter(b,1,demode_signal5_tongbu);
bpf_signal6_tongbu =filter(b,1,demode_signal6_tongbu);
Energy_Detect2 = 0;
 for ii = 1:3500
            Energy_Detect2 = Energy_Detect2 + ((bpf_signal1_tongbu(ii))^2 + (bpf_signal2_tongbu(ii))^2+(bpf_signal3_tongbu(ii))^2+(bpf_signal4_tongbu(ii))^2+(bpf_signal5_tongbu(ii))^2+(bpf_signal6_tongbu(ii))^2)/6;
 end
a = Energy_Detect1/Energy_Detect2;
 if a > 1.01
     tongbu_start = tongbu_start -70;
     tongbu_end = tongbu_start + 6999;
      if tb == 100
%          fprintf('同步失败!\n')
         tongbuchaoshi = tongbuchaoshi + 1;
         Send_Clock = tongbu_end + 1;
         Initial_Send_Clock = Send_Clock;
         buhuo_ok = 0;
         break
     end
 else  if a < 0.99
         tongbu_start = tongbu_start +70;
         tongbu_end = tongbu_start + 6999;
         if tb == 100
%          fprintf('同步失败!\n')
         tongbuchaoshi = tongbuchaoshi + 1;
         Send_Clock = tongbu_end + 1;
         Initial_Send_Clock = Send_Clock;
         buhuo_ok = 0;
         break
         end
     else 
%          fprintf('同步成功!\n')
         tongbu_ok = 1;
         break
      end
 end
end%同步阶段结束

end%大循环结束
% fprintf('大循环结束。\n')
% judge2
% tb
% xujing
% tongbuchaoshi
% Send_Clock
% tongbu_start
%------------测试------------
% tongbu_start = fix((tongbu_start+3500)/7000)*7000
% tongbu_start = 28001;
% judge2 = 5;
%------------跳频解调阶段------------
d = fir1(128,0.0104);%13K低通
e = fir1(512,0.0008);%1K低通
f = fir1(800,[0.003,0.005]);%2.5k-7.5k 0带通
g = fir1(800,[0.007,0.009]);%7.5k-12.5k 1带通

if(buhuo_ok ==1&&tongbu_ok ==1)
receive_signal = FHMsignal(1,tongbu_start:tongbu_start + 1000*N_de-1);%同步后待解跳的100码元信号
% carrier = FHP1;
% signal2 = receive_signal.*carrier;
% signal3 = filter(c,1,signal2);
% signal3W = fft(signal3,length(signal3));
% t = 0:Fs/(length(signal3)-1):Fs;
% plot(t,abs(signal3W));


demode_carrier = [];

% FHP_tiaoxu = [315 455 525 245 385 175 105]*10^3;%Fs/4=625KHz  2FSK第一零点带宽5+2*2.5=10KHz 频谱平移后占25KHz 跳频间隔为70KHz 足够
t2 = 0:1/Fs:(7000-1)/Fs;
FHP1 = cos(2*pi.*FHP_tiaoxu(1).*t2);
FHP2 = cos(2*pi.*FHP_tiaoxu(2).*t2);
FHP3 = cos(2*pi.*FHP_tiaoxu(3).*t2);
FHP4 = cos(2*pi.*FHP_tiaoxu(4).*t2);
FHP5 = cos(2*pi.*FHP_tiaoxu(5).*t2);
FHP6 = cos(2*pi.*FHP_tiaoxu(6).*t2);
FHP7 = cos(2*pi.*FHP_tiaoxu(7).*t2);

% FHP = [315 455 525 245 385 175 105]*10^3;%Fs/4=625KHz  2FSK第一零点带宽5+2*2.5=10KHz 频谱平移后占25KHz 跳频间隔为70KHz 足够
% t2 = 0:1/Fs:(7000-1)/Fs;
% FHP1 = (1/sqrt(2))*cos(2*pi.*FHP(1).*t2)+sin(2*pi.*FHP(1).*t2);
% FHP2 = (1/sqrt(2))*cos(2*pi.*FHP(2).*t2)+sin(2*pi.*FHP(2).*t2);
% FHP3 = (1/sqrt(2))*cos(2*pi.*FHP(3).*t2)+sin(2*pi.*FHP(3).*t2);
% FHP4 = (1/sqrt(2))*cos(2*pi.*FHP(4).*t2)+sin(2*pi.*FHP(4).*t2);
% FHP5 = (1/sqrt(2))*cos(2*pi.*FHP(5).*t2)+sin(2*pi.*FHP(5).*t2);
% FHP6 = (1/sqrt(2))*cos(2*pi.*FHP(6).*t2)+sin(2*pi.*FHP(6).*t2);
% FHP7 = (1/sqrt(2))*cos(2*pi.*FHP(7).*t2)+sin(2*pi.*FHP(7).*t2);



de_carrier_storehouse = [FHP1;FHP2;FHP3;FHP4;FHP5;FHP6;FHP7];%本地载波生成逻辑
N_shang = fix(N_de/7);
N_yu = mod(N_de,7);
for i = 1 : N_shang
    demode_carrier = [demode_carrier de_carrier_storehouse(judge2,:)];
    if judge2 == 7 
        judge2 = 0;
    end
    judge2 = judge2 + 1;
end
if N_yu ~= 0
%     if judge2 == 7 
%         judge2 = 0;
%     end
%     judge2 = judge2 + 1;
demode_carrier = [demode_carrier de_carrier_storehouse(judge2,1:N_yu*1000)];
end

% demode_carrier = [FHP5 FHP6 FHP7 FHP1 FHP2 FHP3 FHP4 FHP5 FHP6 FHP7 FHP1 FHP2 FHP3 FHP4 FHP5(1:2000)];

signal_processmix = receive_signal.*demode_carrier;%混频
signal_process2FSK = filtfilt(d,1,signal_processmix);%低通滤波得2FSK信号
% signal_process2FFT = fft(signal_process2,length(signal_process2));
% t = 0:Fs/(length(signal_process2)-1):Fs;
% plot(t,abs(signal_process2FFT))

t_de = 0:1/Fs:(N_de*1000-1)/Fs;
de_FH0 = cos(2*pi*5*10^3*t_de);
de_FH1 = cos(2*pi*10*10^3*t_de);%2FSK解调本地载波生成
%------------0路解调------------
signal_10 = filtfilt(f,1,signal_process2FSK);%过带通
signal_20 = signal_10.*de_FH0;%混频
signal_30 = filtfilt(e,1,signal_20);%过低通,用于抽样判决
%------------1路解调------------
signal_11 = filtfilt(g,1,signal_process2FSK);%过带通
signal_21 = signal_11.*de_FH1;%混频
signal_31 = filtfilt(e,1,signal_21);%过低通，用于抽样判决
%------------抽样判决------------
judger = abs(signal_31)-abs(signal_30);
result_seq = [];
step_de = 1000;
MAXCLOCK_de = N_de*1000;
CLOCK_de = 0;
while CLOCK_de < MAXCLOCK_de
    CLOCK_de = CLOCK_de + step_de;
    CLOCK_judge = CLOCK_de -500;
    if judger(CLOCK_judge)>0
        result_seq = [result_seq 1];
    else
        result_seq = [result_seq 0];
    end
end
a1 = fix((tongbu_start+500)/1000);%寻找发射端正确序列，tongbu_start偏超500认为无效试验
original_seq = x(a1 + 1:a1 + N_de);%对应发射端比特位
ber_seq = xor(result_seq,original_seq);
be = sum(ber_seq);
c = 0;
if tongbu_ok == 1
    c = abs(mod(tongbu_start+3500,7000)-3500);
    cc = [cc c];
end
% %------------------------------------------------------------------------------------------------------------------------
% figure(5)
% subplot(2,1,1)
% plot(signal_processmix)
% title('混频信号时域')
% W = 0:Fs/(length(signal_processmix)-1):Fs;
% signal_processWmixW = fft(signal_processmix,length(signal_processmix));
% subplot(2,1,2)
% plot(W,abs(signal_processWmixW))
% title('混频信号频域')
% %------------------------------------------------------------------------------------------------------------------------
% figure(6)
% subplot(2,1,1)
% plot(signal_process2FSK)
% title('过低通恢复2FSK信号时域')
% W = 0:Fs/(length(signal_process2FSK)-1):Fs;
% signal_process2FSKW = fft(signal_process2FSK,length(signal_process2FSK));
% subplot(2,1,2)
% plot(W,abs(signal_process2FSKW))
% title('过低通恢复2FSK信号频域')
% %------------------------------------------------------------------------------------------------------------------------
% figure(7)
% subplot(3,2,1)
% plot(signal_10)
% title('0信号过带通时域')
% subplot(3,2,3)
% plot(signal_20)
% title('0信号混频时域')
% subplot(3,2,5)
% plot(signal_30)
% title('0信号过低通时域')
% subplot(3,2,2)
% plot(signal_11)
% title('1信号过带通时域')
% subplot(3,2,4)
% plot(signal_21)
% title('1信号混频时域')
% subplot(3,2,6)
% plot(signal_31)
% title('1信号过低通时域')
% %------------------------------------------------------------------------------------------------------------------------
% figure(8)
% subplot(2,1,1)
% signal_process2FSKW = fft(signal_process2FSK,length(signal_process2FSK));
% W = 0:Fs/(length(signal_process2FSK)-1):Fs;
% plot(W,abs(signal_process2FSKW))
% % plot(receive_signal)
% subplot(2,1,2)
% demode_carrierW = fft(demode_carrier,length(demode_carrier));
% W = 0:Fs/(length(demode_carrier)-1):Fs;
% plot(W,abs(demode_carrierW))
% plot(demode_carrier)
end
end
% er = [er be];
% end
% a = sum(er)/50000