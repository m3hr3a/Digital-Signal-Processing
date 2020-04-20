function [ G, Wflipped ] = flipDTFT( H, W )
N = length(H); 
Wflipped = -W(N:-1:1);
G = H(N:-1:1);
end
