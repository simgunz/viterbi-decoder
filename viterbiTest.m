clc;
clear all;
close all;

% Since the vitdecoder.c mex function keep the full survivors table, matlab crash if
% long generators are used with high input length, because it runs out of memory.
% If the input length is too short for high SNR the avarege BER results zero so the
% graph is truncated.
% Generators [23,35] with mu = 50000 still works with 4GB RAM
% To decode codes corresponding to longer generator, it is possible to use
% vitdecoderFiniteTable.c, that keeps just a finite number of survivors (5*nu).
% The performance in terms of Pbit are worst. The complexity is higher due to
% rough implementation.

%%% TUNABLE PARAMETERS %%%

g1 = 5;   % Generators in octal notation
g2 = 7;

mu = 100000;        % Input length
finiteTable = 0;    % Decode with finite survivors table (5*nu)
enableMex = 1;      % Enable or disable c implementation of the encoder/decoder
iter = 10;          % Number of simulations (10000 for good results)
EbN0step = 0.1;     % Set to 0.5 to speed things up

%%%%%% SIMULATION %%%%%%%%

tr = poly2trellis(length(de2bi(oct2dec(g1))),[g1,g2]);
spec = distspec(tr);

RVit = 1/2;     % Code rate of the convolutional code
RUnc = 1;       % Code rate of uncoded tx

EbN0 = 1:EbN0step:9;
EbN0dB = 10*log10(EbN0);

gammaVit = 2*RVit*EbN0;
gammaUnc = 2*RUnc*EbN0;

if ~finiteTable
    [a,index]=min(abs(qfunc(sqrt(spec.dfree*gammaVit))-5e-6)); % Find the index of the gamma
                                                               % corresponding to Pbit=5e-6
    gammaVit = gammaVit(1:index);   % Compute only for gamma corresponding to
                                    % Pbit > 5e-6 approximatly
end

gammaVitDB = 10*log10(gammaVit);
gammaUncDB = 10*log10(gammaUnc);

PbitVit = zeros(1,length(gammaVit));
PbitUnc = zeros(1,length(gammaUnc));

tic
for k=1:length(gammaVit)
    for i=1:iter
        u_input = round(rand(1,mu));       % Random input sequence
        u_outputVit = convolutionalTxSystem( u_input, g1, g2, enableMex, gammaVitDB(k) ,finiteTable);
        PbitVit(k) = PbitVit(k) + sum(u_input ~= u_outputVit);
    end
    PbitVit(k) = PbitVit(k)/(mu*iter);
end
for k=1:length(gammaUnc)
    for i=1:iter
        u_input = round(rand(1,mu));       % Random input sequence
        u_outputUnc = uncodedTxSystem( u_input, gammaUncDB(k) );
        PbitUnc(k) = PbitUnc(k) + sum(u_input ~= u_outputUnc);
    end
    PbitUnc(k) = PbitUnc(k)/(mu*iter);
end
toc

PbitUncTh = qfunc(sqrt(gammaUnc));
PbitVitTh = qfunc(sqrt(spec.dfree*gammaVit));
% PbitVitTh = bercoding(EbN0,'conv','soft',RVit,spec);

h = figure;
semilogy(EbN0dB,PbitUncTh,'m');
hold;
semilogy(EbN0dB,PbitUnc,'g');
semilogy(EbN0dB(1:length(gammaVit)),PbitVitTh,'c');
semilogy(EbN0dB(1:length(gammaVit)),PbitVit);
line([0,10],[1e-5,1e-5],'Color','r');

legend('Uncoded Thoretical BER','Uncoded Simulated BER','CC Thoretical BER','CC Simulated BER');
xlabel('Eb/N0 [dB]');
ylabel('Pbit');
mkdir('output');
save('output/workspace');
saveas(h,'output/figure');
saveas(h,'output/figure','pdf');