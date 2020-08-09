 function [Iestimated,correlation,similitude,t,dvf] = dvf_demons_estimated (fixedScene,movingScene,pyramLevels,initialItera,smoothLevels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function non rigid registration to estimate the contour
% propagation and displacement vector fields of a tomography sequence, from
% ct and cst structure set, based-demons algorithm
% 
% call
%   [Iestimated,movingScenery,correlation,similitude,t,dvf] = matRad_ContourPropagation(fixedScene,movingScene,structNumber,ct,cst,pyramLevels,initialItera,smoothLevels)
%
% input
%   movingScene:   scene number of the moving image with values between 1
%                  and 10 from ct structure
%   fixedScene:    image scene fixed from ct structure
%
%   imageModality: modality of the image
%   ct:            matRad ct structure 
%   cst:           matRad cst struct 
%   pyramLevels:   number of multi-resolution image pyramid levels to use,
%                  specified as the positive integer scalar
%   initialItera:  number of iterations, specified as a positive integer
%                  scalar 
%   smoothLevels:  smoothing applied at each iteration. This parameter
%                  controls the amount of diffusion-like regularization, it
%                  applies the standard deviation of the Gaussian
%                  smoothing. Values typically are in the range [0.5 , 3.0]
%
% output
%   Iestimated:    estimated image 
%   movingScenery: original scenario
%   correlation:   computes the correlation coefficient
%                  between images in 3D Iestimated and movingScenery image
%                  of cst structure
%   similitude:    computes the SSD metric between images in 3D Iestimated
%                  and movingScenery image of cst structure
%   t:             time in seconds
%   dvf:           displacement vector fields         
%
% References
%   -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    % iterations vector by levels given initial iteration
    N=zeros(pyramLevels);
    N(1) = initialItera;
    for j=2:pyramLevels
        N(j) = round(N(j-1)*0.5);
    end

            % compute displacement vector fields based-demons algorithm, where fixed image is ct.cubeHU{1,fixedScene} and moving image is ct.cubeHU{1,movingScene}
            [dvf,~] = imregdemons(fixedScene,movingScene,N(1:pyramLevels),'PyramidLevels',pyramLevels,'AccumulatedFieldSmoothing',smoothLevels);
            
            Iestimated=imwarp(fixedScene,dvf);
%             movingScenery = movingTomogra;
            correlation = corr3D(movingScene,Iestimated); % correlation coefficient          
            similitude = ssd3D(movingScene,Iestimated);  % metric ssd
            
    t=toc;
end
