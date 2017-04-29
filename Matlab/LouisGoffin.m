    boxL   = 2.
    Lfloor = 0.7
    Lwater = 0.5
    sep = 0.05/4/2 % Rboman = 0.05/4/2
    
    kernel = 'cubic'       
    
    h_0      = 0.2%rboman=0.06/4/2   % initial smoothing length [m]
    c_0      = 35.0       % initial speed of sound  [m/s]
    rho_0    = 1000.0     % initial density [kg/m^3]
    dom_dim  = boxL       % domain size (cube)
    alpha    = 0.5        % artificial viscosity factor 1
    beta     = 0.0        % artificial viscosity factor 2
    T  = 3.0        % simulation time
    saveInt  = 0.01/2     % save interval
    
    % mobile particles
    disp(['Fluid'])
     s=sep
     o=[(boxL-Lwater)/2 (boxL-Lwater)/2 (boxL)/2+0.5]
     L=[Lwater Lwater Lwater]

    
    % fixed particles
    
    % obstacle
     disp(['obstacle'])
     s=sep
     o=[(boxL-Lfloor)/2 (boxL-Lfloor)/2 boxL/2]
     L=[Lfloor Lfloor sep]

    % floor
     disp(['obstacle'])
     s=sep
     o=[0 0 0]
     L=[boxL boxL sep]
    
    % x=0
     disp(['x=0'])
     s=sep
     o=[0 0 2*sep]
     L=[sep boxL boxL-2*sep]

    % y=0
     disp(['y=0'])
     s=sep
     o=[2*sep 0 2*sep]
     L=[boxL-4*sep sep boxL-2*sep]
     
    
    % x=L
     disp(['x=L'])
     s=sep
     o=[boxL-sep 0 2*sep]
     L=[sep boxL boxL-2*sep]

    
    % y=L
     disp(['y=L'])
     s=sep
     o=[2*sep boxL-sep 2*sep]
     L=[boxL-4*sep sep boxL-2*sep]
