%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LoadParameters.m
%
% Returns a column vector of fixed model parameters including:
% - Transition rates (phiA, phiR, phiE)
% - Natural decay rates (deltaA, deltaR, deltaE, deltaP)
% - Proliferation rates (omegaEL, omegaRL)
% - T cell activation rates (lambdaEL, lambdaR)
% - Carrying capacity (C)
% - APC and peptide related parameters (sigma, alpha)
% - Beta cell death rate (kappa)
% - Regulatory T-cell mechanism parameters (r1, r2)
%
% These parameters are used as inputs for the immune cell and beta cell 
% ODE model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function params = LoadParameters()

    % Transition rate parameters
    phiA = 0.9242;
    phiR = 0.05;
    phiE = 0.05;

    % Natural decay parameters
    deltaA = 0.25;
    deltaR = 0.25;
    deltaE = 0.25;
    deltaP = 1;

    % Proliferation Parameters
    omegaEL = 1.735;
    omegaRL = 1.735;

    % T cell activation parameters
    lambdaEL = 1;
    lambdaR = 1;
    
    % Carrying Capacity parameter
    C = 2.5 * 10^7;
    
    % APC and Peptide Parameters
    sigma = 100;
    alpha = 1*10^(-5);
    
    % Beta Death Rate
    kappa = 0.14*10^(-6);

    % Reg T Mechanism Parameters
    r1 = 0;
    r2 = 0;

    % Return parameters as a column vector
    params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaRL C ...
              phiR deltaR sigma alpha deltaP kappa r1 r2]';

end