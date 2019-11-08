function [mutualInfoCrit, mutualInfoMtx] = nearNeiConditionalMutualInformation_for_LA1(dataPreprocessed,idTargets,idConditionalTerm,metric,k,funcDir,homeDir)
% Description:
%
% This function is implemented with MuTE (http://mutetoolbox.guru/) and TSTool
% toolbox (http://www.physik3.gwdg.de/tstool/), to use the nearest-neighbour
% estimator, the original m-file nearNeiFirstTermCondEnt.m in /MuTE/nnnue
% should be substituted. Corresponding functions of the LA method for
% binning and linear estimators may be modified in the same manner, which
% located respectively at /MuTE/commonFunctions/conditionalEntropy.m and
% /MuTE/linue/linearEntropy.m
%


%% Building the delay matrix
mutualInfoMtx                           = buildingEntropyMtx(dataPreprocessed,idTargets,idConditionalTerm);

%% Evaluating the low-dimensional approximation of conditional mutual information for ranking candidate vectors

if size(mutualInfoMtx{1,1},1)>2
    
    y=mutualInfoMtx{1,1}(1,:); %target vector
    V=mutualInfoMtx{1,1}(2:end-1,:); %embedding vectors that already chosen
    x=mutualInfoMtx{1,1}(end,:); %causal driver
    
    Ixy=evalNNMutualInfo([x',y'],metric,k,funcDir,homeDir);
    Ixv=zeros(1,size(V,1));
    Ixvy=zeros(1,size(V,1));
    
    for ii=1:size(V,1)
        Ixv(1,ii)=evalNNMutualInfo([x',V(ii,:)'],metric,k,funcDir,homeDir);
        Ixvy(1,ii)=evalConditionalMutualInfo([x',y',V(ii,:)'],metric,k,funcDir,homeDir);
    end    
    mutualInfoCrit = Ixy - 2*mean(Ixv) + 2*mean(Ixvy);
else
    mutualInfoCrit = evalNNMutualInfo(mutualInfoMtx{1,1}',metric,k,funcDir,homeDir);
end

