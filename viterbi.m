clc;
%clear all;

%%%%%%% TUNABLE PARAMETERS

g1 = [1 0 1];
g2 = [1 1 1];

mu = 200000;                      % Input length
nu = length(g1) - 1;            % Number of memory elements
ns = 2^nu;

Psymbol = 0.05;     % Psymbol = Q(1/std_deviation_w)
                    % Percentage of erroneous transmitted BPAM symbols (symbols received with
                    % opposite sign)

%mexEnabled = 0;     % Enable or disable c implementation of the encoder/decoder
test = 1;
test1 = 0;
iter = 1000;
% Y(u,s,k) Output map, store the vector [y1, y2]
% in position (u,s,:) associated to the input u and the state s
% S(u,s) State update map (matrix), store the state sNext
% in position (u,s) associated to the input u and the state s
% N(s,k,i) Neighbors map, store the k-th neighbor [u, sN] associated to
% the state s in position (s,k,:) 

% For Matlab all the values of input and states are shifted by 1
% To match matlab index convention the range of stored
% state is 1-2^nu and the range of stored input is 1-2
% To match C index convention the range of stored input is 0-1
% and the range of stored state is 0-(2^nu)-1

[Y,S,N] = buildMaps(g1,g2);

if ~mexEnabled       % Translate to MATLAB index system
    S = S + 1;
    N = N + 1;    
end


tic
for i=1:iter
    
    u_input = round(rand(1,mu));       % Input sequence

    if test
        u_input = uTest;        
    end
    if test1
        u_input = [0 1 1 1 0 1 0 0 0];
        mu = 9;
    end
    u = [u_input zeros(1,nu)];         % Input extended with zero padding
    

    if ~mexEnabled
        u = u + 1;                % Match the correct index system
    end
    
    %%%%%%% ENCODER %%%%%%%

    % Output vector, each column is [y1; y2]

    
    if mexEnabled
        y = vitencoder(u,Y,S);
    else
        s = 1;                  % Initial state
        y = zeros(2,mu + nu);
        for k=1:(mu + nu)
            y(:,k) = Y(u(k),s,:);
            s = S(u(k),s);
        end
    end


    %%%%%%%%% P/S %%%%%%%%%

    d = y(:)';


    %%% BPAM MODULATOR %%%%

    sTx = d;
    sTx(sTx==0) = -1;    % Map 0 to -1


    %%%%%%% CHANNEL %%%%%%%

    std_deviation_w = 1/(qfuncinv(Psymbol));            % Noise standard deviation
    r = sTx + randn(1,length(sTx)) * std_deviation_w;   % Received signal
    
     if test
        r = rTest;
     end
    symbolErr = sum(sign(sTx) ~= sign(r));              % Number of erroneous tx symbols


    %%%%%%%%% S/P %%%%%%%%%

    r = reshape(r,2,length(r)/2);


    %%% VITERBI DECODER %%%
    %% ELIMINA RIDIMENSIONAMENTO ARRAY
    if mexEnabled
        u_output = vitdecoder(r,Y,S,N);
    else
        gamma = ones(1,ns)*(-inf);    % Metrics column vector
        gamma(1) = 0;
        gammanew = gamma;
        survivors = zeros(2,mu+nu,ns);       % Survivors array of vectors        
        % M = 
        for i=1:(mu + nu)                   % For every output [y1; y2]
            for j=1:ns                    % For every possible state                
                    
                M(:,1) = squeeze(Y(N(j,1,1),N(j,1,2),:));
                M(:,2) = squeeze(Y(N(j,2,1),N(j,2,2),:));
                M(M==0) = -1;
                tempgamma = r(:,i)'*M + gamma(N(j,1,2):N(j,2,2));
                if tempgamma(1) > tempgamma(2) 
                    bestK = 1;
                else
                    bestK = 2;
                end
                gammanew(j) = tempgamma(bestK);
                survivors(:,i,j) = squeeze(N(j,bestK,:));                
                
            end
            gamma = gammanew - max(gammanew);            
        end
        u_output = zeros(1,mu+nu);
        ss = 1;    
        for k=fliplr(1:(nu+mu))        
            u_output(k) = survivors(1,k,ss);
            ss = survivors(2,k,ss);
        end
        u_output = u_output(1:mu) - 1;      % Eventually translate from matlab index to real value
    end
       
    residualErr = sum(u_input ~= u_output);

end
toc
disp(['Symbols errors = ',num2str(symbolErr)]);
disp(['Residual errors = ',num2str(residualErr)]);

