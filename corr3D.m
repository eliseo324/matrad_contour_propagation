function d = corr3D(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to compute the correlation coefficient between 
% 3D images A and B. Values are in the range [0,1]
% 
% call
%   d = 
%
% input
%   A:       
%   B:              
%
% output                                                                                    
%   d:                         
%
% References
%   -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    A = A(:) - mean2(A(:));
    B = B(:) - mean2(B(:));
    d = sum(sum(A.*B))/sqrt(sum(sum(A.*A))*sum(sum(B.*B)));
       
end