% ditr_complete.m

%%Analyze stability in the complete case for the Domestic (current 
%and expected) Inflation Taylor Rule

clear;
clear all 
global sigma shi chi fi beta eta gama eps sc tita v delta landa k1 k2 w mu si

LOAD_OLD = 0;

if LOAD_OLD == 1
    load ditr_complete
else
    
    %********* Initialize fundamental Parameters **************************************
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
    delta=0; %0.144;    %Openess in production {0,0.144}

    %*********Reduced parameter**********
    tau=sigma*(1-gama)/((1-gama)^2+sigma*gama*eta*(2-gama));
    %New keynesian Phillips
    landa=((1-tita)*(1-tita*beta)/tita)*((1-v)*(1-delta)/(1-v+...
        delta*chi));
    k1=chi+(tau/(1-gama))*((1-v+delta*chi)/((1-v)*(1-delta)));
    %control to reduce to complete market
    k2=0;
    %Is equation
    w=sigma/(sigma-fi);
    mu=((1-gama)/(sigma-fi*(1-gama)))*(1-gama+eta*gama*(2-gama)*sigma/(1-gama));
    si=eta*gama*fi*(2-gama)/((1-gama)*(sigma-fi));

    %******THREE CASES TO COMPARE******
    %we will analyze the imcomplete case, which is the baseline, and then we
    %compare with other two cases: (I) Comlete case (kdos=0), (II) and the
    %close economy (gama=0, and delta=0).

    %**********Stability analysis***********
    decimals=100;    
    fixmax=4;
    fipimax=4;
    fiq=0.6;
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

    for j=1:m+1
        fix=(j-1)/decimals;
        fixs(j)=fix;
        for s=1:n+1
            fipi=(s-1)/decimals;
            fipis(s)=fipi;
            iteration=j*s;
            disp(['Iteration #:', int2str(iteration)])

    %************ Calculate jacobian and eigenvalues***************************
    a11 = (1+fix*(mu-si*(1-gama))+(landa*k1/beta)*(mu-si*(1-gama)))/w;
    a12 = (fipi*(mu-si*(1-gama))-(1/beta)*(mu-si*(1-gama)))/w;
    a21 = -landa*k1/beta;
    a22 = 1/beta;

    A=[a11 a12; a21 a22];
    [eigvec,eigval]=eig(A);
    egv1=eigval(1,1);
    egv2=eigval(2,2);

    eigv1(s,j)=egv1;
    eigv2(s,j)=egv2;

            if abs(real(egv1))>1&&abs(real(egv2))>1                 
               %deter1=[-1,-1]; 
               stable = [fipi, fix];
               stable1 = [stable1; stable];
            else
               deter1=[fipi,fix];           
            end
            deter=[deter;deter1];

        end
    end

    obs1=length(deter);
    M1=sortrows(deter,2); 
    COMB=M1(obs1,:);  %Length of indeterminacy

    save('ditr_complete.mat', 'stable1', 'deter', 'COMB');
end

display('INDETERMINACY REGION FOR DITR')
disp(COMB)

% figure;
% hold on
% 
% subplot(1,1,1)
% plotmatrix(deter(:,1),deter(:,2))
% title('Open economy with complete markets and DITR');
% xlabel('policy reaction to inflation');
% ylabel('policy reaction to output-gap');
% 
% hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%::::::::::::EXPECTED DOMESTIC INFLATION TAYLOR RULE (EDITR):::::::::::::::
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%same notation as before by including "e" at the end

if LOAD_OLD == 1
    load ditre_complete
else
    
    %Creating output matrix
    fixse=zeros(m,1); %parameter of output gap
    fipise=zeros(n,1);   %parameter of inflation gap
    detere=[];  %Give the combinations of fipi and fix that imply determinacy
    stable2 = [];
    
    %eigenvalues
    eigve1=zeros(n,m);
    eigve2=zeros(n,m);

    for j=1:m+1
        fixe=(j-1)/decimals;
        fixse(j)=fixe;
        for s=1:n+1
            fipie=(s-1)/decimals;
            fipise(s)=fipie;
            iteration=j*s;
                disp(['Iteration #:', int2str(iteration)])

    %************ Calculate jacobian and eigenvalues***************************
    ae11 =(w-fixe*(mu-si*(1-gama)))^(-1)*(1+...
        (landa*k1/beta)*(1-fipie)*(mu-si*(1-gama)));
    ae12 =-(w-fixe*(mu-si*(1-gama)))^(-1)*(1/beta)*(1-fipie)*(mu-si*(1-gama));
    ae21 = -landa*k1/beta;
    ae22 = 1/beta;

    AE=[ae11 ae12; ae21 ae22];
    [eigvece,eigvale]=eig(AE);
    egve1=eigvale(1,1);
    egve2=eigvale(2,2);

    eigve1(s,j)=egve1;
    eigve2(s,j)=egve2;

            if abs(real(egve1))>1&&abs(real(egve2))>1               
               %detere1=[-1,-1]; 
               stable = [fipie, fixe];
               stable2 = [stable2; stable]; 
            else
               detere1=[fipie,fixe];         
            end
           detere=[detere;detere1];

        end
    end

    Me1=sortrows(detere,1);
    
save('ditre_complete.mat', 'stable2', 'detere');
end

display('INDETERMINACY REGION FOR FB-DITR')

% figure;
% hold on
% 
% subplot(1,1,1)
% plotmatrix(detere(:,1),detere(:,2))
% title('Open economy with complete markets and FB-DITR');
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

% Patch diagram for DITR case
figure
    title('Open complete-markets economy with DITR');
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
    
    print('-depsc', 'oe-cm-ditr')
    hgsave('oe-cm-ditr')

% Patch diagram for FB-DITR case
figure
    title('Open complete-markets economy with FB-DITR');
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
    print('-depsc', 'oe-cm-fbditr')
    hgsave('oe-cm-fbditr')

