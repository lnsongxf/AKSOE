
% mertrCPI.m

% Analyze stability of the incomplete markets for the CPI (current and expected) 
% managed exchange rate Taylor rule. 


clear all 
global sigma shi chi fi beta eta gama eps sc tita v delta landa k1 k2 w mu si

% USER SETTINGS:
% -------------------------------------------------------------------------
LOAD_OLD = 1; % load saved results

POLICY = 1;   % Set to 1 for MERTR & FBMERTR;  == 0 for CPITR or FBCPITR

if LOAD_OLD == 1
    if POLICY == 1
        load mertr
    elseif POLICY == 0
        load cpitr
    end
else
    
    %********* Initialize fundamental Parameters **************************
    sigma=5;        %IES consumption {1,5}
    shi=-1;         %Share desutility labor
    chi=0.47;       %Frisch elasticity {0.47,3}
    fi=10^(-6);  %Endogenerity of discounting {0,10^(-6)}
    beta=0.99;      %stationary discounting
    raro=0;         %parameter to calibrate discounting
    eta=1.5;        %Elasticity substitution in composite consumption {1,1.5}
    gama=0.4;       %share in composite consumption (openess in demand) {0,0.4}
    eps=6;          %Elasticity substitution in agregator consumption
    sc=0.8;         %Consumption to output ratio
    tita=0.75;      %Calvo probability
    v=-2;           %Elasticity of substititution in production {0,-2}
    delta=0.144;    %Openess in production {0,0.144}

    %*********Reduced parameter**********
    %New keynesian Phillips
    landa=((1-tita)*(1-tita*beta)/tita)*((1-v)*(1-delta)/(1-v+...
        delta*chi));
    k1=(sigma/(1-gama))+chi;
    %control to reduce to complete market
    k2=delta*(1-v+chi)/((1-gama)*(1-v)*(1-delta))+gama/(1-gama)-...
        sigma*eta*gama*(2-gama)/((1-gama)^2);
    %Is equation
    w=sigma/(sigma-fi);
    mu=((1-gama)/(sigma-fi))*(1-gama+eta*gama*(2-gama)*(sigma-fi)/(1-gama));
    si=eta*gama*fi*(2-gama)/((1-gama)*(sigma-fi));

    %******THREE CASES TO COMPARE******
    %we will analyze the imcomplete case, which is the baseline, and then we
    %compare with other two cases: (I) Comlete case (kdos=0), (II) and the
    %close economy (gama=0, and delta=0).

    %**********Stability analysis***********
    decimals=101;    
    fixmax=4;
    fipimax=4;
    
    % Reaction to ER
    if POLICY == 1
        fis=0.61;   % MERTR / FBMERTR
    elseif POLICY == 0
        fis = 0;    % CPITR / FBCPITR
    end

    
    m=fixmax*decimals;  %interaciones en fix 
    n=fipimax*decimals; %interaciones en fipi 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %::::::::::::DOMESTIC INFLATION TAYLOR RULE (DITR)::::::::::::::::::::::::
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Creating the output matrix
    fixs=zeros(m,1); %parameter of output gap
    fipis=zeros(n,1);   %parameter of inflation gap
    deter=[];  %Give the combinations of fipi and fix that imply determinacy
    stable1 = [];
    
    %eigenvalues
    eigv1=zeros(n,m);
    eigv2=zeros(n,m);
    eigv3=zeros(n,m);
    eigv4=zeros(n,m);
    count=zeros(n,m);
    for j=1:m+1
        fix=(j-1)/decimals;
        fixs(j)=fix;
        for s=1:n+1
            fipi=(s-1)/decimals;
            fipis(s)=fipi;
                iteration=j*s;
            disp(['Iteration #:', int2str(iteration)])

    %************ Calculate jacobian and eigenvalues***************************
    fx=fix;
    fpi=fipi+fis;
    fs=(gama*fipi+fis)/(1-gama);

    B=[w mu si 0; 0 beta 0 0; 0 (1-gama) 1 0; 0 0 0 1];
    d11=1+mu*fx;
    d12=mu*fpi;
    d13=mu*fs;
    d14=-mu*fs;
    d21=-landa*k1;
    d22=1;
    d23=-landa*k2;
    d24=0;
    d31=(1-gama)*fx;
    d32=(1-gama)*fpi;
    d33=1+(1-gama)*fs;
    d34=-(1-gama)*fs;
    d41=0;
    d42=0;
    d43=1;
    d44=0;
    D=[d11 d12 d13 d14; d21 d22 d23 d24; d31 d32 d33 d34; d41 d42 d43 d44];

    f=B\D;  %This matrix coincides with matrix A
    [eigvec,eigval]=eig(f);
    egv1=eigval(1,1);
    egv2=eigval(2,2);
    egv3=eigval(3,3);
    egv4=eigval(4,4);
    eigv1(s,j)=egv1;
    eigv2(s,j)=egv2;
    eigv3(s,j)=egv3;
    eigv3(s,j)=egv4;
    counter=0;
    %Let's count the number of unstable eigenvalues
        if abs(real(egv1))>1
            counter=counter+1;
        end  
        if abs(real(egv2))>1
            counter=counter+1;
        end
        if abs(real(egv3))>1
            counter=counter+1;
        end  
        if abs(real(egv4))>1
            counter=counter+1;
        end  

        count(s,j)=counter;
    %Let's analyze stability
            if counter==3                
               stable = [fipi, fix];
               stable1 = [stable1; stable];
            else
               deter1=[fipi,fix];   
               deter=[deter; deter1];           
            end

        end
    end

    obs1=length(deter);
    M1=sortrows(deter,2); 
    COMB=M1(obs1,:);  %Length of indeterminacy
    
    if POLICY == 1
        save('mertr.mat', 'stable1', 'deter', 'COMB');
    elseif POLICY == 0
        save('cpitr.mat', 'stable1', 'deter', 'COMB');
    end
end

display('INDETERMINACY REGION FOR MERTR')
disp(COMB)

% figure;
% hold on
% 
% subplot(1,1,1)
% plotmatrix(deter(:,1),deter(:,2))
% title('Open economy with incomplete markets and MERTR');
% xlabel('policy reaction to inflation');
% ylabel('policy reaction to output-gap');
% 
% hold off


if LOAD_OLD == 1
    if POLICY == 1
        load mertre
    elseif POLICY == 0
        load cpitre
    end
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %::::::: Forecast-based Managed exchange rate Taylor Rule (FB-MERTR):
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %same notation as before by including "e" at the end

    %Creating output matrix
    fixse=zeros(m,1); %parameter of output gap
    fipise=zeros(n,1);   %parameter of inflation gap
    detere=[];  %Give the combinations of fipi and fix that imply determinacy
    stable2 = [];

    %eigenvalues
    eigve1=zeros(n,m);
    eigve2=zeros(n,m);
    eigve3=zeros(n,m);
    for j=1:m+1
        fixe=(j-1)/decimals;
        fixse(j)=fixe;
        for s=1:n+1
            fipie=(s-1)/decimals;
            fipise(s)=fipie;
            iteration=j*s;
            disp(['Iteration #:', int2str(iteration)])

    %************ Calculate jacobian and eigenvalues***************************
    fxe=fixe;
    fpie=fipie+fis;
    fse=(gama*fipie+fis)/(1-gama);

    be11=w-mu*fxe;
    be12=mu*(1-fpie);
    be13=si-mu*fse;
    be21=0;
    be22=beta;
    be23=0;
    be31=-(1-gama)*fxe;
    be32=(1-gama)*(1-fpie);
    be33=1-(1-gama)*fse;
    BE=[be11 be12 be13; be21 be22 be23; be31 be32 be33];

    de11=1; 
    de12=0;
    de13=-mu*fse;
    de21=-landa*k1;
    de22=1;
    de23=-landa*k2;
    de31=0;
    de32=0;
    de33=1-(1-gama)*fse;
    DE=[de11 de12 de13; de21 de22 de23; de31 de32 de33];

    fe=BE\DE;
    [eigvece,eigvale]=eig(fe);    
    egve1=eigvale(1,1);
    egve2=eigvale(2,2);
    egve3=eigvale(3,3);
    eigv1(s,j)=egve1;
    eigv2(s,j)=egve2;
    eigv3(s,j)=egve3;

            if abs(real(egve1))>1&&abs(real(egve2))>1&&abs(real(egve3))>1                 
               stable = [fipie, fixe];
               stable2 = [stable2; stable];
            else
               detere1=[fipie,fixe];   
               detere=[detere; detere1];         
            end
        end
    end

    Me1=sortrows(detere,2);

    if POLICY == 1
        save('mertre.mat', 'stable2', 'detere', 'COMB');
    elseif POLICY == 0
        save('cpitre.mat', 'stable2', 'detere', 'COMB');
    end
end

% display('INDETERMINACY REGION FOR FB-DITR')
% figure;
% hold on
% 
% subplot(1,1,1)
% plotmatrix(detere(:,1),detere(:,2))
% title('Open economy with incomplete markets and FB-DITR');
% xlabel('policy reaction to inflation');
% ylabel('policy reaction to output-gap');
% 
% hold off

%==========================================================================
%     GENERATE PATCH SURFACES
%==========================================================================

colortheme1 = [.69  .69  .69];
colortheme2 = [.078  .745  .945];

xmax = 4;
ymax = 4;

% Patch diagram
figure
if POLICY == 1
    title('Open incomplete-markets economy with MERTR');
elseif POLICY == 0
    title('Open incomplete-markets economy with CPITR');
end
    xlabel('\phi_{\pi}');
    ylabel('\phi_{x}');
    axis([0 xmax 0 ymax])
    
    x = stable1(:,1);
    y = stable1(:,2);
    k = convhull(x,y);
    hold on
    
        rectangle('Position',[0,0,xmax,ymax],...
                    'FaceColor',colortheme1);
                
        p1 = patch(x(k), y(k),'w');
        
        set(p1, 'FaceColor',    colortheme2, ...
                'EdgeColor',    colortheme2, ...
                'LineWidth',    1 );
                
    hold off;
    
    if POLICY == 1
        filename = 'oe-im-mertr';
    elseif POLICY == 0
        filename = 'oe-im-cpitr';
    end
    
    print('-depsc', filename)
    hgsave(filename)

% Patch diagram for FB-DITR case
figure
if POLICY == 1
    title('Open incomplete-markets economy with FB-MERTR');
elseif POLICY == 0
    title('Open incomplete-markets economy with FB-CPITR');
end
    xlabel('\phi_{\pi}');
    ylabel('\phi_{x}');
    axis([0 xmax 0 ymax])
    
    
    x = stable2(:,1);
    y = stable2(:,2);
    k = convhull(x,y);
    
    hold on
            
        rectangle('Position',[0,0,xmax,ymax],...
                    'FaceColor',colortheme1);
        
        p2 = patch(x(k), y(k),'w');
        
        set(p2, 'FaceColor',    colortheme2, ...
                'EdgeColor',    colortheme2, ...
                'LineWidth',    1 );
                
    hold off;
    
    % Save this figure as .EPS and .FIG
    if POLICY == 1
        filename = 'oe-im-fbmertr';
    elseif POLICY == 0
        filename = 'oe-im-fbcpitr';
    end

    print('-depsc', filename)
    hgsave(filename)
