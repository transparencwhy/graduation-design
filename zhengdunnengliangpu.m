%����ɲ�����������ƿ��Ӹ������������ͼ
load FHMsignal&FHMsignalW;
Fs = 2.5*10^6;%����Ƶ��
Rs = 2.5*10^3;%��Ԫ��������
N = 560;%50��Ԫ��1��Ԫ1000�����
N_de = 100;%������Ԫ��
Len = N*Fs/Rs;%�ܵ�������
x = randi([0 1],1,N);%���������56����Ϊ������Ԫ
t = 0:1/Fs:(Len-1)/Fs;
menxian = 29;
FHP_local = [313 453 523 243 383 173 103]*10^3;%�����ز��뷢����ز������Ƶ2Khz
b = fir1(30,0.0104);
gaosibaizaosheng = wgn(1,Len,10);
f_total = [];
% function f = Energy_Detect(signal_cut,signal_move)
st = 5001;
ed = st + 3499;
t1 = (0:3499)/Fs;
signal_move = cos(2*pi*FHP_local(3)*t1);
signal_cut = FHMsignal;
% i = 0;
for i = 1:49000
demode_signal1 = signal_cut(1,st-2+i:ed-2+i).*signal_move; 
demode_signal2 = signal_cut(1,st-1+i:ed-1+i).*signal_move; 
demode_signal3 = signal_cut(1,st+i:ed+i).*signal_move; 
demode_signal4 = signal_cut(1,st+1+i:ed+1+i).*signal_move; 
demode_signal5 = signal_cut(1,st+2+i:ed+2+i).*signal_move; 
demode_signal6 = signal_cut(1,st+3+i:ed+3+i).*signal_move; 

bpf_signal1 =filter(b,1,demode_signal1);
bpf_signal2 =filter(b,1,demode_signal2);
bpf_signal3=filter(b,1,demode_signal3);
bpf_signal4 =filter(b,1,demode_signal4);
bpf_signal5 =filter(b,1,demode_signal5);
bpf_signal6 =filter(b,1,demode_signal6);
f = 0;
%  for ii = 1:1000
%             f = f +(bpf_signal3(ii))^2;
%  end
 for ii = 1:3500
            f = f + [(bpf_signal1(ii))^2 + (bpf_signal2(ii))^2+(bpf_signal3(ii))^2+(bpf_signal4(ii))^2+(bpf_signal5(ii))^2+(bpf_signal6(ii))^2]/6;
 end
f_total = [f_total f];
end
plot(f_total(1:49000))
title('����ģ�����������')
ylabel('����')
xlabel('���������ʼ��/��')
