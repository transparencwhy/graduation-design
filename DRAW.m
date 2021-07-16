%画图

BPSKGAILV_pu_you = [99.35 99.33 99.22 99.25 99.22];
JSR = [0 0.5 1 1.5 2];

BPSKGAILV_pu_wu = [99.35 98.95 96.24 70.71 40.31];
JSR = [0 0.5 1 1.5 2];

figure(1);
plot(JSR,BPSKGAILV_pu_wu,'--or',JSR,BPSKGAILV_pu_you,'^g:');
legend('未规避干扰','规避干扰');
title('部分频带阻塞干扰下的基于谱估计的捕获方法捕获概率图');
ylabel('捕获概率/%');
xlabel('干信比');
axis([0 2 0 100]);

BPSKGAILV_chuan_you = [87.35 87.05 86.83 87.06 86.27];
JSR = [0 0.5 1 1.5 2];

BPSKGAILV_chuan_wu = [87.35 84.76 73.9 66.05 52.8 47.99 41.45 33.56 32.79 33.13 33.24];
JSR2 = [0 0.5 0.55 0.58 0.59 0.6 0.65 0.75 1 1.5 2];

figure(2)
plot(JSR2,BPSKGAILV_chuan_wu,'--or',JSR,BPSKGAILV_chuan_you,'^g:');
legend('未规避干扰','规避干扰');
title('部分频带阻塞干扰下的快速扫描式捕获方法捕获概率图');
ylabel('捕获概率/%');
xlabel('干信比');
axis([0 2 0 100]);
%%单音干扰s3
DANYINGAILV_pu_you = [99.35 99.32 99.19 99.23 99.24];
JSR = [0 0.5 1 1.5 2];

DANYINGAILV_pu_wu = [99.35  99.17 76.19 49.01 21.41 14.60 9.48 0 0 0 0];
JSR1 = [0 0.02 0.03 0.04 0.05 0.06 0.07 0.5 1 1.5 2];

figure(3)
plot(JSR1,DANYINGAILV_pu_wu,'--or',JSR,DANYINGAILV_pu_you,'^g:');
legend('未规避干扰','规避干扰');
title('单音干扰下的基于谱估计的捕获方法捕获概率图');
ylabel('捕获概率/%');
xlabel('干信比');
axis([0 2 0 100]);

DANYINGAILV_chuan_you = [86.66 86.97 86.83 86.14 86.66 86.15 84.79];
JSR1 = [0.25 0.5 1 1.5 2 3 5];

DANYINGAILV_chuan_wu = [85.24 73.81 68.93 58.14 50.64 50.30 49.61 49.18 48.60 48.16];
JSR2 = [0.25 0.27 0.28 0.29 0.5 1 1.5 2 3 5]; 
figure(4)
plot(JSR2,DANYINGAILV_chuan_wu,'--or',JSR1,DANYINGAILV_chuan_you,'^g:');
legend('未规避干扰','规避干扰');
title('单音干扰下的快速搜索扫描式捕获方法捕获概率图');
ylabel('捕获概率/%');
xlabel('干信比');
axis([0.25 5 0 100]);