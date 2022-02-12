clear;

load("ber_1_16qam.mat");

load("ber_2_16qam.mat");

% Removal of data outliers
for i = 1:size(ber1,1)
    for j = 1:size(ber1,2)
        if(ber1(i,j) > 0.4)
            ber1(i,j) = 0;
        end
        if(ber2(i,j) > 0.4)
            ber2(i,j) = 0;
        end
    end
end

var1 = zeros(1, size(ber1,1));
sum1 = zeros(1, size(ber1,1));

var2 = zeros(1, size(ber1,1));
sum2 = zeros(1, size(ber1,1));

for k = 1:size(ber1,1)    
        var1(k) = var(ber1(k,:));
        var2(k) = var(ber2(k,:));
    
        sum1(k) = sum(ber1(k,:));
        sum2(k) = sum(ber2(k,:));
end

figure(1)
plot([10:20],var1);
title('Variance of BER with increasing SNR for 16QAM, 1x4 case');
xlabel('SNR in dB'); ylabel('Variance');

figure(2)
plot([10:20],var2);
title('Variance of BER with increasing SNR for 16QAM, 2x2 case');
xlabel('SNR in dB'); ylabel('Variance');

figure(3)
plot([10:20],sum1/size(ber1,1));
title('Average of BER with increasing SNR for 16QAM, 1x4 case');
xlabel('SNR in dB'); ylabel('Average');

figure(4)
plot([10:20],sum2/size(ber1,1));
title('Average of BER with increasing SNR for 16QAM, 2x2 case');
xlabel('SNR in dB'); ylabel('Average');
