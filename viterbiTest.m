clc;
clear all;
close all;

%%%%%%% TUNABLE PARAMETERS

g1 = [1 0 1];
g2 = [1 1 1];

mu = 1000;                      % Input length

Psymbol = 0.05;     % Psymbol = Q(1/std_deviation_w)
                    % Percentage of erroneous transmitted BPAM symbols (symbols received with
                    % opposite sign)

mex = 1;     % Enable or disable c implementation of the encoder/decoder
iter = 10;



gammaDB = 0:0.1:12.5;

PbitVit = zeros(1,length(gammaDB));
PbitUnc = zeros(1,length(gammaDB));    

tic
for k=1:length(gammaDB)    
    for i=1:iter

        u_input = round(rand(1,mu));       % Input sequence

        u_outputVit = convolutionalTxSystem( u_input, g1, g2, mex, gammaDB(k) );
        u_outputUnc = uncodedTxSystem( u_input, gammaDB(k) );

        PbitVit(k) = PbitVit(k) + sum(u_input ~= u_outputVit);
        PbitUnc(k) = PbitUnc(k) + sum(u_input ~= u_outputUnc); 
    end
    PbitVit(k) = PbitVit(k)/(mu*iter);
    PbitUnc(k) = PbitUnc(k)/(mu*iter);
end
toc

PbitUncTh = qfunc(sqrt(10.^(gammaDB/10)));
PbitVitTh = qfunc(sqrt(5*10.^(gammaDB/10)));
figure;
semilogy(10*log10((10.^(gammaDB/10))/2),PbitUncTh,'r');
hold;
semilogy(10*log10((10.^(gammaDB/10))/2),PbitUnc,'g');
semilogy(10*log10((10.^(gammaDB/10))/5),PbitVitTh,'b');
semilogy(10*log10((10.^(gammaDB/10))/5),PbitVit);
