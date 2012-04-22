function [ symbol ] = M( y )
%M Summary of this function goes here
%   Detailed explanation goes here

    y(find(y==0))=-1;
    symbol = y;
end

