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

    

params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaRL C ...
          phiR deltaR sigma alpha deltaP kappa r1 r2]';

end