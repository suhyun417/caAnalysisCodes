function [C, Cparams] = computeSlidingWindowCorrAcrossTrials(matsdf, sizeW)

nMovie = length(matsdf); %size(matsdf,3);
C=struct;

for iMovie = 1:nMovie
    
    tempMatSDF = matsdf{iMovie}; %matsdf(:,:,iMovie);
    
    % parameters for sliding window
    % sizeW = 5000; % msec
    stepW = sizeW/5;
    numW = (size(tempMatSDF,1)-sizeW)./stepW+1; %(size(matsdf,1)-sizeW)./stepW+1;
    % edStW = size(S(iC).matsdf,1)-sizeW;
    
    setPair = combnk(1:size(tempMatSDF,2), 2); % set of all possible pairs of trials
    numPair = size(setPair,1);
    
    
    for iW = 1:numW
        for iP=1:numPair
            [R,P,RLO,RUP]=corrcoef(tempMatSDF(stepW*(iW-1)+1:stepW*(iW-1)+sizeW, setPair(iP,:)));
            C(iMovie).R(iP,iW) = R(1,2);
            C(iMovie).P(iP,iW) = P(1,2);
            C(iMovie).RLO(iP,iW) = RLO(1,2);
            C(iMovie).RUP(iP,iW) = RUP(1,2);
        end
    end
    
end

Cparams.sizeW = sizeW;
Cparams.stepW = stepW;
Cparams.numW = numW;
Cparams.setPair = setPair;
Cparams.numPair = numPair;

        
        
