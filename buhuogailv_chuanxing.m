%���в��������֤����������Ӧ���ϣ�������������������²��������Ҫ49000+7000+7000*3=77000�㼴77����Ԫ����77�����˸�105���ɣ����ݵ���105*1000=105000��
clc
clear all
% load BPSKforchuanxing
load BPSKforpu;
% load FHMsignal&FHMsignalW;
er = [];
gailv = [];
for iiii = 1:100
%------------��ʼ��������------------
Fs = 2.5*10^6;%����Ƶ��
Rs = 2.5*10^3;%��Ԫ��������
N = 105;%50��Ԫ��1��Ԫ1000�����
N_de = 50;%������Ԫ��
Len = N*Fs/Rs;%�ܵ�������
x = randi([0 1],1,N);%���������56����Ϊ������Ԫ
t = 0:1/Fs:(Len-1)/Fs;
menxian = 35;
%------------�����ź���ɢ��------------
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
% % % % % figure(1)
% % % % % % subplot(1,1,1)
% % % % % plot(t1,signal);
% % % % % axis([0 50 -2 2]);xlabel('��Ԫ');ylabel('����');title('Դ��Ϣ����')
%------------2FSK����------------
f0 = 5*10^3;
f1 = 10*10^3;%2FSK���ز�Ƶ��f0=5KHz,f1=10KHz
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
% % % % % figure(2)
% % % % % subplot(2,1,1);
% % % % % plot(t,FSK);title('2FSK�����ź�ʱ����');xlabel('ʱ��(��)');ylabel('����');
% % % % % axis([0 0.002 -2 2]);
% % % % % FSKW = fft(FSK,Len);
% % % % % t1 = 0:Fs/(Len-1):Fs;
% % % % % subplot(2,1,2);
% % % % % plot(t1,abs(FSKW))
% % % % % title('2FSKƵ��ͼ');xlabel('Ƶ��(Hz)');ylabel('����');
% % % % % axis([0 15000 -inf +inf]);

%------------Ƶ�ʵ��ʼ��------------
FHPlocal = [313 453 523 243 383 173 103]*10^3;%�����ز��뷢����ز������Ƶ2Khz
FHP_tiaoxu = [315 455 525 245 385 175 105]*10^3;%Fs/4=625KHz  2FSK��һ������5+2*2.5=10KHz Ƶ��ƽ�ƺ�ռ25KHz ��Ƶ���Ϊ70KHz �㹻
FHP_shunxu = [105 175 245 315 385 455 525]*10^3;%Fs/4=625KHz  2FSK��һ������5+2*2.5=10KHz Ƶ��ƽ�ƺ�ռ25KHz ��Ƶ���Ϊ70KHz �㹻
% % a1 = 200;b1 = 350;
% % a2 = 400; b2 = 550;
% % a3 = 600; b3 = 750;
% % a4 = 800; b4 = 950;
% % a5 = 1000; b5 = 1150;
% % a6 = 1200; b6 = 1350;
% % a7 = 1400; b7 = 1550;
% a_seq = [200 400 600 800 1000 1200 1400];
% b_seq = [350 550 750 950 1150 1350 1550];
% % aa1 = 1488; bb1 = 1641;
% % aa2 = 1600; bb2 = 1750;
% % aa3 = 1712; bb3 = 1862;
% % new_aa = [aa1 aa2 aa3];
% % new_bb = [bb1 bb2 bb3];
% aa_seq = [1488 1600 1712];
% bb_seq = [1638 1750 1862];
% 
% % c1 = 500;d1 = 670;
% % c2 = 900;d2 = 1070;
% % c3 = 1300; d3 = 1470;
% % c4 = 1700; d4 = 1870;
% % c5 = 2100; d5 = 2270;
% % c6 = 2500; d6 = 2670;
% % c7 = 2900; d7 = 3070;
% c_seq = [500 900 1300 1700 2100 2500 2900];
% d_seq = [670 1070 1470 1870 2270 2670 3070];

% % cc1 = 3076; dd1 = 3246;
% % cc2 = 3300; dd2 = 3470;
% % cc3 = 3524; dd3 = 3694;
% % new_cc = [cc1 cc2 cc3];
% % new_dd = [dd1 dd2 dd3];
% cc_seq = [3076 3300 3524];
% dd_seq = [3246 3470 3694];
%------------����Ƶ����������------------
noise_filter = fir1(512,[0.08,0.16]);%100K��200KHz
%------------��˹����������------------
% gaosibaizaosheng = wgn(1,Len,10);
%------------������Ų���------------
t_noise = 0:1/Fs:(Len-1)/Fs;
s1 = sin(2*pi*FHP_shunxu(1)*t_noise);
s2 = sin(2*pi*FHP_shunxu(2)*t_noise);
s3 = sin(2*pi*FHP_shunxu(3)*t_noise);
s4 = sin(2*pi*FHP_shunxu(4)*t_noise);
s5 = sin(2*pi*FHP_shunxu(5)*t_noise);
s6 = sin(2*pi*FHP_shunxu(6)*t_noise);
s7 = sin(2*pi*FHP_shunxu(7)*t_noise);
%------------����������ѡ��------------
% interference_or_noise = s1;
ganrao =  zeros(1,105000) ;
%------------�ŵ����------------
FFTW = 0:Fs/(10000-1):Fs;
ganrao_clip = ganrao(1,1:10000);
ganrao_clipW = fft(ganrao_clip,10000);
% plot(FFTW,abs(ganrao_clipW));
ganrao_clipW = abs(ganrao_clipW);
volume_jiance = [max(ganrao_clipW(400:450)) max(ganrao_clipW(680:730)) max(ganrao_clipW(960:1010)) max(ganrao_clipW(1240:1290)) max(ganrao_clipW(1520:1570)) max(ganrao_clipW(1800:1850)) max(ganrao_clipW(2080:2130))];
volume_jiance = volume_jiance>20000000;
volume_beiyong = [595 665 735 805 875 945 1015]*10^3;%�߸�����Ƶ��
volume_beiyong_zhongpin = [593 663 733 803 873 943 1013]*10^3;
for i = 1:7%����Ӧ��ܸ���Ƶ�㡣
    if volume_jiance(i) == 1
        b = find(FHP_shunxu == (35+70*i)*10^3);
        FHP_shunxu(b) = volume_beiyong(1);
        b = find(FHP_tiaoxu == (35+70*i)*10^3);
        FHP_tiaoxu(b) = volume_beiyong(1);
        FHPlocal(b) = volume_beiyong_zhongpin(1);
        volume_beiyong(1:6) = volume_beiyong(2:7);%����
        volume_beiyong_zhongpin(1:6) = volume_beiyong_zhongpin(2:7);%����
%         a_seq(i) = aa_seq(1);
%         aa_seq(1:2) = aa_seq(2:3);%����
%         b_seq(i) = bb_seq(1);
%         bb_seq(1:2) = bb_seq(2:3);%����
%         c_seq(i) = cc_seq(1);
%         cc_seq(1:2) = cc_seq(2:3);%����
%         d_seq(i) = dd_seq(1);
%         dd_seq(1:2) = dd_seq(2:3);%����
    end
end

%------------m��������------------
% an = [1 0 0];%3λ�Ĵ���
% len = length(an);
% L = 2^len-1;%m���г���Ϊ7
% m = zeros(1,L)
% m(1)= an(3);
% for i = 2:L
%     an0=mod((an(1)+an(3)),2)
%     an(2:len)=an(1:len-1);%����1λ
%     an(1)=an0;
%     m(i) = an(3);
% end
%------------��Ƶ�ز���Ӧ------------
% FHP = [105 175 245 315 385 455 525]*10^3;%Fs/4=625KHz  2FSK��һ������5+2*2.5=10KHz Ƶ��ƽ�ƺ�ռ25KHz ��Ƶ���Ϊ70KHz �㹻
t2 = 0:1/Fs:(7000-1)/Fs;
FHP1 = cos(2*pi.*FHP_shunxu(1).*t2);
FHP2 = cos(2*pi.*FHP_shunxu(2).*t2);
FHP3 = cos(2*pi.*FHP_shunxu(3).*t2);
FHP4 = cos(2*pi.*FHP_shunxu(4).*t2);
FHP5 = cos(2*pi.*FHP_shunxu(5).*t2);
FHP6 = cos(2*pi.*FHP_shunxu(6).*t2);
FHP7 = cos(2*pi.*FHP_shunxu(7).*t2);

% FHP = [205 275 345 415 485 555 625]*10^3;%Fs/4=625KHz  2FSK��һ������5+2*2.5=10KHz Ƶ��ƽ�ƺ�ռ25KHz ��Ƶ���Ϊ70KHz �㹻
% t2 = 0:1/Fs:(7000-1)/Fs;
% FHP1 =  (1/sqrt(2))*(cos(2*pi.*FHP(1).*t2)+sin(2*pi.*FHP(1).*t2));
% FHP2 =  (1/sqrt(2))*(cos(2*pi.*FHP(2).*t2)+sin(2*pi.*FHP(2).*t2));
% FHP3 = (1/sqrt(2))*(cos(2*pi.*FHP(3).*t2)+sin(2*pi.*FHP(3).*t2));
% FHP4 = (1/sqrt(2))*(cos(2*pi.*FHP(4).*t2)+sin(2*pi.*FHP(4).*t2));
% FHP5 = (1/sqrt(2))*(cos(2*pi.*FHP(5).*t2)+sin(2*pi.*FHP(5).*t2));
% FHP6 = (1/sqrt(2))*(cos(2*pi.*FHP(6).*t2)+sin(2*pi.*FHP(6).*t2));
% FHP7 = (1/sqrt(2))*(cos(2*pi.*FHP(7).*t2)+sin(2*pi.*FHP(7).*t2));
%------------�ź���Ƶ����------------
carrier = [];
MAXCLOCK = N*Fs/Rs;%����ʱ��t����ѡȡ��Ϊ50000
CLOCK = 0;
step = 7000; %��56����Ԫ��ÿ��Ԫ0.0004s��7����Ԫһ��Ƶ��1/0.0004/7��357��ÿ��
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
    an(2:len)=an(1:len-1);%����1λ
    an(1)=an0;
end
%------------��Ӹ���------------
FHMsignal = awgn(FHMsignal,-10);
%------------������Ƶͼ��------------ 
%------------------------------------------------------------------------------------------------------------------------
% % % % % figure(3)
% % % % % plot(pattern/10^3,'s','markerfacecolor','r','markersize',18); %'s'���飬'markerfacecolor'������,'r'����ɫ��'markersize'���
% % % % % %��С
% % % % % ylabel('kHz');
% % % % % xlabel('��Ԫ');
% % % % % title('��Ƶͼ��');
% % % % % axis([1 56 0 700]);
% % % % % grid on;
%------------��Ƶ�ź�ʱ��Ƶ�����------------
% yes = 0;
% if FSK.*carrier_total == FHMsignal
%     yes = 1;
% else 
%     yes = 0;
% end
%------------------------------------------------------------------------------------------------------------------------
% % % % % figure(4)
% % % % % FHMsignalfftW = fft(FHMsignal,length(FHMsignal));
% % % % % t1 = 0:Fs/(Len-1):Fs;
% % % % % subplot(2,1,1);plot(t,FHMsignal);xlabel('ʱ��(s)');ylabel('����');title('��Ƶ�����ź�ʱ����');
% % % % % axis([0 0.0032 -inf +inf]);
% % % % % subplot(2,1,2);plot(t1/10^3,abs(FHMsignalfftW));xlabel('Ƶ�ʣ�KHz��');ylabel('����');title('��Ƶ�����ź�Ƶ����');

%------------��ʱ����Ҷ�任������Ƶ�ź�------------
% [S W T] = mf_spectrogram(FHMsignal,7000,6000,2.5*10^6);
% my_pcolor( W , T , S );

% spectrogram(FHMsignal,1000,500,1024,2.5*10^6,'yaxis')%�����źţ����β����������������β����ص�������fft����������Ƶ�ʣ��任������
%------------��������������-------------
% Fs = 2.5*10^6;%����Ƶ��
% Rs = 2.5*10^3;%��Ԫ��������
% N = 560;%50��Ԫ��1��Ԫ1000�����
% Len = N*Fs/Rs;%�ܵ�������
% x = randi([0 1],1,N);%���������56����Ϊ������Ԫ
% t = 0:1/Fs:(Len-1)/Fs;
% t2 = 0:Fs/(1000-1):Fs;
%------------������������------------
% FHMsignal = awgn(FHMsignal,1);
% y = wgn(1,Len,10);
% FHMsignal = FHMsignal + y;
% FHMsignal = FHMsignal + cos(2*pi*205*10^3*t);
%------------�±�Ƶ����Ƶ׼��------------
t1 = 0:1/Fs:(1000-1)/Fs;
% FHPlocal = [103 243 523 453 383 173 313]*10^3;%�����ز��뷢����ز������Ƶ2Khz
% % local_carrier = cos(2*pi*FHPlocal(2).*t1);
% % adddmax = 1000;
% % ii = 0;
% % for i = 1 : 500
% local_carrierW = fft(local_carrier,length(local_carrier));
% plot(t2,abs(local_carrierW));
%------------��Ƶ�ź��뱾���ز��źŻ�Ƶ����Ƶ�ź�------------
% % demode_signal = FHMsignal(7001-500+i:7001-500+i+999).*local_carrier;
% %------------�˲�ǰʱ��------------
% figure(1);
% subplot(2,1,1);plot(demode_signal);title('�˲�ǰʱ��');
% %------------�˲�ǰƵ��------------
% demode_signalW = fft(demode_signal,length(demode_signal));
% subplot(2,1,2);plot(t2,abs(demode_signalW));title('�˲�ǰƵ��');
% axis([0 10^5 -inf +inf]);
% %------------ǰ1000����Ƶ�ź�Ƶ��------------
% figure(2);
% FHMsignalW = fft(FHMsignal(1:1000),1000);
% plot(t2,abs(FHMsignalW));title('ǰ1000����Ƶ�ź�Ƶ��');
%------------��Ƶ�źŹ��˲���------------
b = fir1(48,0.0104);
% % bpf_signal = filter(b,1,demode_signal);
% %------------�˲���ʱ��------------
% figure(3);
% subplot(2,1,1);plot(bpf_signal);title('�˲���ʱ��');
% %------------�˲���Ƶ��------------
% bpf_signalW = fft(bpf_signal,length(bpf_signal));
% subplot(2,1,2);plot(t2,abs(bpf_signalW));title('�˲���Ƶ��');
%------------ƽ���ۼ��������------------
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
%------------����Ԥ�趨------------
Send_Clock = randi([1,49000],1,1);
Initial_Send_Clock = Send_Clock;
MAX_CLOCK = Len; 
Send_STEP = 7000;
judge1 = 0;
license = 0;
T_buhuo = 0;
buhuo_ok = 0;
loubu = 0;
queren_ok = 0;
tongbu_ok = 0;
xujing = 0;
tongbuchaoshi = 0;
time = 0;
judge1 = randi([1,7]);
%------------�����߼�------------
while ( buhuo_ok == 0)%�����ʼ����.
    Send_Clock = Send_Clock + Send_STEP;
%     TEMP_signal1 = FHMsignal(1,(Send_Clock - Send_STEP + 1):(Send_Clock));
        while (1)
            time = time + 1;
        if time == 8
            time = 0;
            loubu = 1;%������©��
%             fprintf('����©��!\n')
            break
        end
        if(judge1 == 7)
            judge1 = 0;
        end
        judge1 = judge1 + 1;
        TEMP_STEP = time * 1000;
%         TEMP_signal2 = TEMP_signal1(1,(TEMP_STEP - 1000 +1):(TEMP_STEP));
local_carrier_buhuo = cos(2*pi*FHPlocal(judge1).*t1);
demode_signal1_buhuo = FHMsignal(1,(Send_Clock - Send_STEP + 1)+(TEMP_STEP - 1000) - 2 :(Send_Clock - Send_STEP + 1)+(TEMP_STEP-1)- 2).*local_carrier_buhuo; 
demode_signal2_buhuo = FHMsignal(1,(Send_Clock - Send_STEP + 1)+(TEMP_STEP - 1000) - 1:(Send_Clock - Send_STEP + 1)+(TEMP_STEP-1)  - 1).*local_carrier_buhuo; 
demode_signal3_buhuo = FHMsignal(1,(Send_Clock - Send_STEP + 1)+(TEMP_STEP - 1000):(Send_Clock - Send_STEP + 1)+(TEMP_STEP-1)).*local_carrier_buhuo; 
demode_signal4_buhuo = FHMsignal(1,(Send_Clock - Send_STEP + 1)+(TEMP_STEP - 1000) + 1:(Send_Clock - Send_STEP + 1)+(TEMP_STEP-1) + 1).*local_carrier_buhuo; 
demode_signal5_buhuo = FHMsignal(1,(Send_Clock - Send_STEP + 1)+(TEMP_STEP - 1000) + 2:(Send_Clock - Send_STEP + 1)+(TEMP_STEP-1) + 2).*local_carrier_buhuo; 
demode_signal6_buhuo = FHMsignal(1,(Send_Clock - Send_STEP + 1)+(TEMP_STEP - 1000) + 3:(Send_Clock - Send_STEP + 1)+(TEMP_STEP-1) + 3).*local_carrier_buhuo; 

% b = fir1(48,[0.0008,0.0024]);
bpf_signal1_buhuo =filter(b,1,demode_signal1_buhuo);
bpf_signal2_buhuo =filter(b,1,demode_signal2_buhuo);
bpf_signal3_buhuo=filter(b,1,demode_signal3_buhuo);
bpf_signal4_buhuo =filter(b,1,demode_signal4_buhuo);
bpf_signal5_buhuo =filter(b,1,demode_signal5_buhuo);
bpf_signal6_buhuo =filter(b,1,demode_signal6_buhuo);
Energy_Detect = 0;
 for ii = 1:1000
            Energy_Detect = Energy_Detect + ((bpf_signal1_buhuo(ii))^2 + (bpf_signal2_buhuo(ii))^2+(bpf_signal3_buhuo(ii))^2+(bpf_signal4_buhuo(ii))^2+(bpf_signal5_buhuo(ii))^2+(bpf_signal6_buhuo(ii))^2)/6;
 end

%         mix_signal = TEMP_signal2.*local_carrier;
%         bpf_signal = filter(b,1,mix_signal);
%         Energy_Detect = 0;
%             for i = 1:1000
%             Energy_Detect = Energy_Detect + (bpf_signal(i))^2;
%             end
%             
            if Energy_Detect > menxian;%��ʼ��������
            buhuo_ok = 1;
            if (fix(((Send_Clock - Send_STEP + 1)+(TEMP_STEP - 1000))/7000)==judge1-1||fix(((Send_Clock - Send_STEP + 1)+(TEMP_STEP - 1000))/7000)==judge1-1+7)
                license = 1;
            end
            T_buhuo = mod((Send_Clock-Initial_Send_Clock) + 1,7000);%1Ϊ����������1�������һ�γ�ʼ������̷�����©��������ʾ©�����س̶�   
%             fprintf('��ʼ����!\n')
            break
            end

        end
        
        if buhuo_ok == 1
            break
        end
        if loubu ==1
            break
        end
end%�˳���ʼ����
    
while(buhuo_ok == 1)%���Խ��벶��ȷ��
judge2 = judge1;
%     ST_local = Send_Clock - Send_STEP + 1 +(judge - 1) * 1000 - 500-3500;
%     ED_local = Send_Clock - Send_STEP + 1 +(judge - 1) * 1000 - 500 + 3500 - 1;%ȷ�����ش���ʼ����λ
for i = 1 : 3
    check_start = Send_Clock - Send_STEP + 1 +(time - 1 ) * 1000 + i* 7000;
    check_end = Send_Clock - Send_STEP + 1 +(time) * 1000-1 + i * 7000;

%     check_signal = FHMsignal(1,check_start:check_end);
    if judge2 == 7
        judge2 = 0;
    end
    judge2 = judge2 + 1;
    local_carrier_buhuoqueren = cos(2*pi*FHPlocal(judge2).*t1);
    
    
demode_signal1_buhuoqueren = FHMsignal(1,check_start-2+i:check_end-2+i).*local_carrier_buhuoqueren; 
demode_signal2_buhuoqueren = FHMsignal(1,check_start-1+i:check_end-1+i).*local_carrier_buhuoqueren; 
demode_signal3_buhuoqueren = FHMsignal(1,check_start+i:check_end+i).*local_carrier_buhuoqueren; 
demode_signal4_buhuoqueren = FHMsignal(1,check_start+1+i:check_end+1+i).*local_carrier_buhuoqueren; 
demode_signal5_buhuoqueren = FHMsignal(1,check_start+2+i:check_end+2+i).*local_carrier_buhuoqueren; 
demode_signal6_buhuoqueren = FHMsignal(1,check_start+3+i:check_end+3+i).*local_carrier_buhuoqueren; 

bpf_signal1_buhuoqueren =filter(b,1,demode_signal1_buhuoqueren);
bpf_signal2_buhuoqueren =filter(b,1,demode_signal2_buhuoqueren);
bpf_signal3_buhuoqueren=filter(b,1,demode_signal3_buhuoqueren);
bpf_signal4_buhuoqueren =filter(b,1,demode_signal4_buhuoqueren);
bpf_signal5_buhuoqueren =filter(b,1,demode_signal5_buhuoqueren);
bpf_signal6_buhuoqueren=filter(b,1,demode_signal6_buhuoqueren);
Energy_Detect = 0;
 for ii = 1:1000
            Energy_Detect = Energy_Detect + ((bpf_signal1_buhuoqueren(ii))^2 + (bpf_signal2_buhuoqueren(ii))^2+(bpf_signal3_buhuoqueren(ii))^2+(bpf_signal4_buhuoqueren(ii))^2+(bpf_signal5_buhuoqueren(ii))^2+(bpf_signal6_buhuoqueren(ii))^2)/6;
 end
 
 
 
    if Energy_Detect < menxian%�����ظ�ȷ������
        break
    end
    if i == 3
        queren_ok = 1;
    end
end

if queren_ok == 1
%     fprintf('ȷ�ϳɹ�!\n')
    break
else
%     fprintf('�����龯!\n')
    xujing = xujing + 1;
    buhuo_ok = 0;
    Send_Clock = check_end + 1;
    Initial_Send_Clock = Send_Clock;
    break
end
             
end%�����ظ�ȷ�Ͻ׶ν���
if (queren_ok == 1&&license ==1)
    gailv = [gailv 1];
else
    gailv = [gailv 0];
end
end
summ = sum(gailv);
summ/10000
% 
% T_buhuo
% judge1
% judge2
% tongbu_start
% tb
% xujing
% loubu
% tongbuchaoshi
