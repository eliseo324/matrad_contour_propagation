function d = ssd3D(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to compute the similitude coefficient between 
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

%     C=A+B;
%     intercepSize = length(find(A + B == 2));
%     intercepSize = 2*sum(abs(C(:)));
    diff = imabsdiff(A, B);
    d = sum(sum(diff(:).^2));
    
%      vol_diff = (imdef - imfix);
%         err = sum(sum(sum(sum(sum(vol_diff.^2))))) * dvol/2;
%         g = vol_diff * dvol;
%         perr = vol_diff * dvol;
       
end