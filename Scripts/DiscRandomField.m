function RF = DiscRandomField(COORD,RFinput)
%  
%      RF = DiscRandomField(COORD, RFinput)
%      Compute data for random field discretization
%  
%  
%fprintf('\n * Discretizing the random field ...');

switch RFinput.Type
 case 'Gaussian'
  RF = DiscGaussianRandomField(COORD,RFinput);
  
 case 'Lognormal'
  %%------------------------------------------------------------%%
  %      compute the mean and Stdv of the underlying Gaussian field
  %
  %      NB : RF.Mean and RF.Stdv contain the values pertaining
  %           to the underlying GAUSSIAN field
  %
  %           RF.LNMean and RF.LNStdv contain the real mean value
  %           and standard deviation of the lognormal field
  %%------------------------------------------------------------%%
  RFinput.Stdv = sqrt(log((RFinput.LNStdv / RFinput.LNMean)^2 +1));
%   RFinput.Mean =  log(RFinput.LNMean) - 0.5 * (RFinput.Stdv)^2;
  RFinput.Mean =  log(RFinput.LNMean^2/sqrt(RFinput.LNMean^2 + RFinput.LNStdv^2));
  
  RF = DiscGaussianRandomField(COORD,RFinput);
  RF.LNMean = RFinput.LNMean;
  RF.LNStdv = RFinput.LNStdv;
  
 otherwise
end
