% MQERTR_OPTIMAL.M
%--------------------------------------------------------------------------

% Analyze stability of the incomplete markets model under the family of  
% forward looking rules of the form:
%
% i(t) = fi_pi*E(pi(t+1)) + fi_x*E(x(t+1)) + fi_q*E(q(t+1))
%
% This class of rules contains the optimal policy under a Markovian 
% time-consistent optimal rule. 
%
% NEED TO CHECK OPTIMAL RULE POINT!
%
% NOTE: This is the limit incomplete market economy with 
%       (i  ) delta = 0 (or alpha = 1)
%       (ii ) nu --> 0
%       (iii) fi --> 0
% So that we have the tractable derivation of second order welfare
% consistent loss function.


clear all 
global sigma shi chi fi beta eta gama eps ...
                                sc tita v delta landa k1 k2 w mu si

% USER SETTINGS:
% -------------------------------------------------------------------------
LOAD_OLD = 0; % load saved results

WAIT_BAR = 0;
    
    % Creating output matrix: Pre-allocate memory space

        %decimals = 26; 

        fixmin = -5;
        fixmax = 5;

        fipimin = -5;
        fipimax = 5;

        m = 101;
        n = 101;
        
    %********* Initialize fundamental Parameters **************************
    sigma=5;    % IES consumption {1,5}
    shi=-1;     % Share disutility labor
    chi=0.47;   % Frisch elasticity {0.47,3}
    fi=10^(-6); % Endogenerity of discounting {0,10^(-6)}
    beta=0.99;  % stationary discounting
    raro=0;     % parameter to calibrate discounting
    eta=1.5;    % Elasticity substitution in composite consumption {1,1.5}
    gama=0.4;   % share in composite consumption (openess in demand){0,0.4}
    eps=6;      % Elasticity substitution in agregator consumption
    sc=0.8;     % Consumption to output ratio
    tita=0.75;  % Calvo probability
    v= -2;       % Elasticity of substititution in production {0,-2}
    delta=0.144;    %0.144;% Openess in production {0,0.144}
    
    %******** Model composite parameters **********************************
        
        % New keynesian Phillips
        k = ((1-tita)*(1-tita*beta)/tita);
        
        landa = k*((1-v)*(1-delta) /(1-v+delta*chi));
        k1=(sigma/(1-gama))+chi;

        % control to reduce to complete market
        k2=delta*(1-v+chi)/((1-gama)*(1-v)*(1-delta))+gama/(1-gama)...
                                    - sigma*eta*gama*(2-gama)/((1-gama)^2);

        % IS equation
        w = sigma/(sigma-fi);
        mu = ((1-gama)/(sigma-fi))*(1-gama+eta*gama*(2-gama)...
                                                    *(sigma-fi)/(1-gama));
        si=eta*gama*fi*(2-gama)/((1-gama)*(sigma-fi));
    
    %******* Calculate optimal time-consistent FB-MERTR representation ****
    
    % Quadratic loss function weights: (R_y, R_e, R_pi) ...
    
    iota = eps/(eps-1); % Steady state markup
    
    R_pi = eps/(sc*iota*k) + eps*(1+chi)*(1+chi*eps)/(2*k); % on (pi_{H})^2
    
    w_denominator = sc*mu*(gama + sigma*(eta-1)  ...
                       + chi*(gama + sc*eta*(1-sigma) + sigma*(eta-1)) );
    
    w_1 = (gama + sigma*eta*(1-sc)+(iota-1)*(1-sc*eta)) / w_denominator;
    
    w_2 = (sigma + (iota-1) + iota*chi) / w_denominator;
    
    w_3 = 0;
    
    w_4 = ( gama + (iota-1) + eta*(sigma+sc*chi*iota) ) / w_denominator;
    
    A_y = [ chi*(2+chi)      sigma       -1      0   ;
            sigma           -sigma^2    sigma    0   ;
            -1              sigma        -1      0   ;
            0                0            0      0   ];
    
    D_y = [ 0                0                      0              0   ;
            0     sc*(1-gama)*(1-(1-gama)*sc)   -eta*sc*(1-sc)     0   ;
            0              -eta*sc*(1-sc)       (eta^2)*sc*(1-sc)  0   ;
            0                0                      0              0   ];
    
    G_y = zeros(4,4);         G_y(4,4) = -gama*(1-eta)*(1-gama*(1-gama));
    
    H_y = [ -1               sigma                     -1          0   ;
            sigma    -sigma*((1-sc)*sigma+1+sc)         sigma      0   ;
            -1               sigma                     -1          0   ;
            0                0                      0              0   ];
    
    M_y = zeros(4,4);  M_y(1,1) = (1+chi)/(iota*sc); M_y(2,2) = -(1-sigma);
    
    R_y = M_y + w_1*A_y + w_2*D_y + w_3*G_y + (w_4/sc)*H_y;
    
    % Target state to co-state affine mapping ...
    
    N_y = [ 1,                               0           ;
            1/((1-gama)*sc),      -eta*gama/((1-gama)^2) ;
            0,                    -gama/(1-gama)         ;
            0,                               1           ];
        
    % Loss function in terms of target states x(t) ...
    
    S_y = (N_y') * R_y * N_y;
    
        sy_11 = S_y(1,1);       % weight on y(t)^2
        sy_12 = S_y(1,2);       % weight on y(t) * q(t)
        sy_21 = S_y(2,1);       % weight on q(t) * y(t)
        sy_22 = S_y(2,2);       % weight on q(t)^2
    
    % FOC (trade-off) for optimal time consistent policy 
    f_pi = landa*R_pi*(mu*k1/(1-gama) -k2);
    f_x = sy_12 + mu*sy_11/(1-gama);
    f_q = sy_22 + mu*sy_21/(1-gama);  
    
    % Composite parameters in optimal rule
    p_pi = beta - landa*( k1*mu - k2*(1-gama) );
    p_x = landa*k1*w;
    p_q = landa*(k1*si + k2);
    p_i = landa*( k1*mu + k2*(1-gama) );
    
    opt_denominator = (f_q*(1-gama) + f_x*mu + f_pi*p_i);
    
    % Optimal time consistent \phi_{pi} ...
    fipi_opt = ( f_q*(1-gama) - f_x*mu + f_pi*p_pi ) / opt_denominator;
    
    % Optimal time consistent \phi_{x} ...
    fix_opt = ( f_x*w + f_pi*p_x ) / opt_denominator;
    
    % Optimal time consistent \phi_{q} ...
    fis_opt = ( f_q + f_x*si + f_pi*p_q ) / opt_denominator;
    
    disp('Optimal Markovian policy, (fi*) = (fipi,fix,fis) is')
    disp([fipi_opt,fix_opt,fis_opt])
    
    
    
    % Set up mesh for (fipi, fix)
    % Control this \phi_q := fis ...
     
    fis_SET = sort([fis_opt, 0.6, 0.8, 0.9, 1.5, 2]);
    
    %fis_SET = sort([fis_SET, fis_opt]);
    
        Ncoord = length(fis_SET);     

        fipise = linspace(fipimin,fipimax,n);  % parameter of inflation gap
        fixse = linspace(fixmin,fixmax,m);     % parameter of output gap

        detere = inf*ones(m*n, 2);        % (fipi,fix) => Indeterminacy
        stable2 = inf*ones(m*n, 2);       % (fipi,fix) => Determinacy
        unstable = inf*ones(m*n, 2);      % (fipi,fix) => Nonexistent REE
        
        %eigenvalues
        eigve1 = -inf*ones(n,m);
        eigve2 = -inf*ones(n,m);
        eigve3 = -inf*ones(n,m);

if LOAD_OLD == 0   
    
    for choose_fis = 1:Ncoord
        
            fis = fis_SET(choose_fis); % Fix the response to E{q(t+1)}

        % MODEL
        %------------------------------------------------------------------

        

        % MERTR-optimal: basically like FB-MERTR without lagged q(t-1) term
        % same notation as before by including "e" at the end


        % Generate combinations in set of (fipi,fix):
        fispace = allcomb(fipise, fixse); % Needs ALLCOMB.M
            LP = length(fispace);
        
        if WAIT_BAR == 1  
            h = waitbar(0,['Constructing Case #',int2str(choose_fis),...
                                                       ' Please wait...']);
        end
        
        for j = 1 : LP % Loop over pairs (fipi,fix)

            fipie = fispace(j,1);
            fixe = fispace(j,2);
            
            %************ Calculate jacobian and eigenvalues***************
            fxe = fixe;
            fpie = fipie;
            fse = fis; %(gama*fipie+fis)/(1-gama);

            be11 = w - mu*fxe;
            be12 = mu*(1 - fpie);
            be13 = si-mu*fse;
            be21 = 0;
            be22 = beta;
            be23 = 0;
            be31 = -(1-gama)*fxe;
            be32 = (1-gama)*(1-fpie);
            be33 =1-(1-gama)*fse;

            BE = [  be11 be12 be13; 
                    be21 be22 be23; 
                    be31 be32 be33  ];

            de11 = 1; 
            de12 = 0;
            de13 = 0; % -mu*fse; No lag term q(t) on optimal rule
            de21 = -landa*k1;
            de22 = 1;
            de23 = -landa*k2;
            de31 = 0;
            de32 = 0;
            de33 = 0; %1-(1-gama)*fse;

            DE = [  de11 de12 de13; 
                    de21 de22 de23; 
                    de31 de32 de33  ];

            FE = BE\DE;

            [eigvece,eigvale] = eig(FE);  

            egve1=eigvale(1,1);
            egve2=eigvale(2,2);
            egve3=eigvale(3,3);

            if abs(real(egve1)) > 1 && abs(real(egve2)) > 1 ...
                                            && abs(real(egve3)) > 1                 
               stable = [fipie, fixe];
               stable2(j,:) = stable;
            elseif abs(real(egve1)) <= 1 || abs(real(egve2)) <= 1 ...
                                            || abs(real(egve3)) <= 1
               detere1=[fipie,fixe];   
               detere(j,:) = detere1;
            elseif abs(real(egve1)) < 1 && abs(real(egve2)) < 1 ...
                                            && abs(real(egve3)) < 1
               unstable1=[fipie,fixe];   
               unstable(j,:) = unstable1;
            end
            percent = j/LP;
        end
        
      
        
        stable2 = stable2( stable2 ~= Inf );
        stable2 = reshape(stable2, length(stable2)/2, 2);

        detere = detere( detere ~= Inf );
        detere = reshape(detere, length(detere)/2, 2);
        
        unstable = unstable( unstable ~= Inf );
        unstable = reshape(unstable, length(unstable)/2, 2);

        save(strcat('mqertr_optimal_',int2str(choose_fis),'.mat'), ...
                                           'fis_SET','stable2', 'detere');
        
        if WAIT_BAR == 0
            disp(['Case #', int2str(choose_fis), ': DONE!',...
                    ' phi_{q} = ',num2str(fis_SET(choose_fis))])
        end
    end
    
    if WAIT_BAR == 1
       close(h)
    end
end


%==========================================================================
%     GENERATE PATCH SURFACES
%==========================================================================

figure

    %title('Open IM economy with MERTR-OPTIMAL type');
    
    
    colortheme1 = [1  1  1];             % Percentage RGB: white
    
    
    colortheme2 = { [.929  .929  .929], ...
                    [.671  .671  .671], ...
                    [.571  .571  .571], ...
                    [.471  .471  .471], ...
                    [.4    .4      .4], ...
                    [.3    .3      .3]      };  
                                         % Percentage RGB: grayscale
                                         
    if length(colortheme2) < Ncoord
        disp('ERROR: Please specify more RGB vectors for colortheme2')
        break
    end

    xmin = fipimin;
    xmax = fipimax;

    ymin = fixmin;
    ymax = fixmax;
            
    axis([xmin xmax ymin ymax])
    
    rectangle('Position',[xmin,ymin,xmax-xmin,ymax-ymin],...
                                'FaceColor',colortheme1);
    plot(fipi_opt,fix_opt,'or','MarkerFaceColor','r')
    
    opt_str = char(sprintf('(%5.2g,%5.2g,%5.2g)',...
                                    fipi_opt,fix_opt,fis_opt));
    
    %text(fipi_opt,1.05*fix_opt,sprintf('(%5.2g,%5.2g,%5.2g)',...
    %                                        fipi_opt,fix_opt,fis_opt))
    gtext(['(\phi^{\ast}) = ',opt_str]) 
    
    if WAIT_BAR == 1
        h1 = waitbar(0,'Constructing Plots. Please wait...');
    end
    
    for choose_fis = 1:Ncoord % Loop over 'fis' settings ...

        if choose_fis == 1
            load mqertr_optimal_1;     
        elseif choose_fis == 2
            load mqertr_optimal_2;  
        elseif choose_fis == 3
            load mqertr_optimal_3;
        elseif choose_fis == 4
            load mqertr_optimal_4;
        elseif choose_fis == 5
            load mqertr_optimal_5;
        end
                
        if isempty(stable2) ~= 1 || isempty(unstable) ~= 1
                
                % Then do Patch diagram ...

                if isempty(stable2) ~= 1
                    % Stable/Unique REE region
                    x = stable2(:,1);
                    y = stable2(:,2);
                    k = convhull(x,y);
                    
                    hold on

                    p1 = patch(x(k), y(k),'w');

                    set(p1, 'FaceColor',    colortheme2{choose_fis}, ...
                            'EdgeColor',    colortheme2{choose_fis}, ...
                            'LineWidth',    1 );
                    pause(3)
                    hold off;
                end
                
                if isempty(unstable) ~= 1
                    % Unstable: nonexistent REE region
                    x1 = stable2(:,1);
                    y1 = stable2(:,2);
                    k1 = convhull(x1,y1);
                    
                    hold on

                    p2 = patch(x1(k), y1(k),'g');

                    set(p2, 'FaceColor',    colortheme2{choose_fis}, ...
                            'EdgeColor',    colortheme2{choose_fis}, ...
                            'LineWidth',    1 );
                    pause(3)
                    hold off;
                    
                end
        else
            disp(['Case phi_s= ',num2str(fis_SET(choose_fis)),...
                         ' WARNING: All economies have indeterminate REE'])     
        end
        xlabel('\phi_{\pi}');
        ylabel('\phi_{x}');
        
        % Generate auto EPS graphics files
                filename = 'oe-im-mqertr-opt';
                print('-depsc', filename)
                hgsave(filename)
        
        if WAIT_BAR == 1
            waitbar(choose_fis/Ncoord,h1)
        end
    end
    
    if WAIT_BAR == 1
        close(h1)
    end