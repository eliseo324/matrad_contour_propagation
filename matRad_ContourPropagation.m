 function [Iestimated,movingScenery,correlation,similitude,t,dvf] = matRad_ContourPropagation(fixedScene,movingScene,structNumber,ct,cst,pyramLevels,initialItera,smoothLevels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function non rigid registration to estimate the contour
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
%   strucNumber:   structure number cst; 2 for complete tomography   
%                  3 for liver struct, 4 for PTV struct
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
    
    fixedTomogra = ct.cubeHU{1,fixedScene}; %fixed tomography 
    movingTomogra = ct.cubeHU{1,movingScene}; %moving tomography 
    
    % compute displacement vector fields based-demons algorithm, where fixed image is ct.cubeHU{1,fixedScene} and moving image is ct.cubeHU{1,movingScene}
    [dvf,~] = imregdemons(fixedTomogra,movingTomogra,N(1:pyramLevels),'PyramidLevels',pyramLevels,'AccumulatedFieldSmoothing',smoothLevels);

    switch structNumber
        case 2
            Iestimated=imwarp(fixedTomogra,dvf);
            movingScenery = movingTomogra;
            correlation = corr3D(movingScenery,Iestimated); % correlation coefficient          
            similitude = ssd3D(movingScenery,Iestimated);  % metric ssd
        case {3,4}  
            % fixed scene 
            cube_fixedScena = zeros(ct.cubeDim);
            structFixed_cst = cst{structNumber,4}{fixedScene};
            [x2,y2,z2] = ind2sub(ct.cubeDim,structFixed_cst);
            for i=1:length(x2)
                cube_fixedScena(x2(i),y2(i),z2(i)) = structFixed_cst(i); 
            end

            % moving scene 
            movingScenery = zeros(ct.cubeDim);
            structMoving_cst = cst{structNumber,4}{movingScene};
            [x1,y1,z1] = ind2sub(ct.cubeDim,structMoving_cst);
            for i=1:length(x1)
                movingScenery(x1(i),y1(i),z1(i)) = structMoving_cst(i);
            end

            % binary case image estimated with the imwarp() function
            Iestimated = imwarp(cube_fixedScena,dvf);
            correlation = corr3D(movingScenery,Iestimated); % correlation coefficient          
            similitude = ssd3D(movingScenery,Iestimated);  % metric ssd
           
    end 
       
    t=toc;
end
