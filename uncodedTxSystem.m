function [ u_output ] = uncodedTxSystem( u, gammaDB )
%CONVOLUTIONALTXSYSTEM Summary of this function goes here
%   Detailed explanation goes here

%%% BPAM MODULATOR %%%%

sTx = u;
sTx(sTx==0) = -1;    % Map 0 to -1


%%%%%%% CHANNEL %%%%%%%

r = awgn(sTx,gammaDB);              % Received signal


%%%%%%% DECODER %%%%%%%

u_output = sign(r);
u_output(u_output==-1) = 0;
    
end

