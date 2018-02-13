% ditr.m

%%Analyze stability our model for the Domestic (current and expected) 
%Inflation Taylor Rule and for the exchange rate Taylor rule. We
%have to manage for the value of fiq:
%- DITR: fiq=0
%- ERTR: fiq>0

clear;
clear all 
global sigma shi chi fi beta eta gama eps sc tita v delta landa k1 k2 ...
                                                                    w mu si

LOAD_OLD = 0; % load saved results ('diter.mat' and 'ditre.mat')

if LOAD_OLD == 1
    load ditr
else
    
    
    %********* Initialize fundamental Parameters **************************
    sigma=1;        %IES consumption {1,5}
    shi=-1;         %Share desutility labor
    chi=0.47;       %Frisch elasticity {0.47,3}
    fi=10^(-6);     %Endogenerity of discounting {0,10^(-6)}
    beta=0.99;      %stationary discounting
    raro=0;         %parameter to calibrate discounting
    eta=1.5;        %Elasticity substitution in composite consumption {1,1.5}
    gama=0.4;       %share in composite consumption (openess in demand) {0,0.4}
    eps=6;          %Elasticity substitution in agregator consumption
    sc=0.8;         %Consumption to output ratio
    tita=0.75;      %Calvo probability
    v=-2;           %Elasticity of substititution in production {0,-2}
    delta=0%; 0.144;    %Openess in production {0,0.144}

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
    decimals=100;    
    fixmax=4;
    fipimax=4;
    fiq=0.6;
    m=fixmax*decimals;  %interaciones en fix 
    n=fipimax*decimals; %interaciones en fipi 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %	      DOMESTIC INFLATION TAYLOR RULE (DITR)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Creating the output matrix
    fixs=zeros(m,1);  %parameter of output gap
    fipis=zeros(n,1); %parameter of inflation gap
    deter=[];         %Give the combinations of fipi and fix that imply determinacy
    stable1 = [];

    %eigenvalues
    eigv1=zeros(n,m);
    eigv2=zeros(n,m);
    eigv3=zeros(n,m);
    for j=1:m+1
        fix=(j-1)/decimals;
        fixs(j)=fix;
        for s=1:n+1
            fipi=(s-1)/decimals;
            fipis(s)=fipi;
            iteration=j*s;
            disp(['Iteration #:', int2str(iteration)])

        %************ Calculate jacobian and eigenvalues***********************
        a11 = (1+fix*(mu-si*(1-gama))+(landa*k1/beta)*(mu-si*(1-gama)))/w;
        a12 = (fipi*(mu-si*(1-gama))-(1/beta)*(mu-si*(1-gama)))/w;
        a13 = (-si+(fiq+(landa*k2/beta))*(mu-si*(1-gama)))/w;
        a21 = -landa*k1/beta;
        a22 = 1/beta;
        a23 = -landa*k2/beta;
        a31 = (1-gama)*(fix-a21);
        a32 = (1-gama)*(fipi-a22);
        a33 = 1+(1-gama)*(fiq-a23);

        A=[a11 a12 a13; a21 a22 a23; a31 a32 a33];
        [eigvec,eigval]=eig(A);
        egv1=eigval(1,1);
        egv2=eigval(2,2);
        egv3=eigval(3,3);
        eigv1(s,j)=egv1;
        eigv2(s,j)=egv2;
        eigv3(s,j)=egv3;
            if abs(real(egv1))>1&&abs(real(egv2))>1&&abs(real(egv3))>1
               stable = [fipi, fix];
               stable1 = [stable1; stable];
               %deter1=[-1,-1]; 
            else
               deter1=[fipi,fix];
               deter=[deter;deter1];
            end


        end
    end


    
    obs1=length(deter);
    M1=sortrows(deter,2); 
    COMB=M1(obs1,:);  %Length of indeterminacy
    
    % SAVE if LOAD_OLD == 0
    save('ditr.mat', 'stable1', 'deter', 'COMB');
    
    display('INDETERMINACY REGION FOR DITR')
    disp(COMB)

end % if LOAD_OLD

    
% figure;
%     hold on
% 
%     subplot(1,1,1)
%     plotmatrix(deter(:,1),deter(:,2))
%     title('Open economy with incomplete markets and DITR');
%     xlabel('\phi_{\pi}');
%     ylabel('\phi_{x}');
% 
%     hold off


if LOAD_OLD == 1
        load ditre
else
    %======================================================================
    %       EXPECTED DOMESTIC INFLATION TAYLOR RULE (EDITR) 
    %======================================================================
    %same notation as before by including "e" at the end

    %Creating output matrix
    fixse = zeros(m,1);     %parameter of output gap
    fipise = zeros(n,1);    %parameter of inflation gap
    detere = [];            %Give the combinations of fipi and fix that 
                            % imply determinacy
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

    %************ Calculate jacobian and eigenvalues***********************
    ae11 =(w-fixe*(mu-si*(1-gama)))^(-1)*(1+...
                            (landa*k1/beta)*(1-fipie)*(mu-si*(1-gama)));
    ae12 =-(w-fixe*(mu-si*(1-gama)))^(-1)*(1/beta)...
                                *(1-fipie)*(mu-si*(1-gama));
    ae13 =(w-fixe*(mu-si*(1-gama)))^(-1)*(-si+...
        (landa*k2/beta)*(1-fipie)*(mu-si*(1-gama)));
    ae21 = -landa*k1/beta;
    ae22 = 1/beta;
    ae23 = -landa*k2/beta;
    ae31 = (1-gama)*(fixe*ae11-ae21*(1-fipie));
    ae32 = (1-gama)*(fixe*ae12-ae22*(1-fipie));
    ae33 = 1+(1-gama)*(fixe*ae13-ae23*(1-fipie));

    AE=[ae11 ae12 ae13; ae21 ae22 ae23; ae31 ae32 ae33];
    [eigvece,eigvale]=eig(AE);
    egve1=eigvale(1,1);
    egve2=eigvale(2,2);
    egve3=eigvale(3,3);
    eigve1(s,j)=egve1;
    eigve2(s,j)=egve2;
    eigve3(s,j)=egve3;
            if abs(real(egve1))>1&&abs(real(egve2))>1&&abs(real(egve3))>1                 
               %detere1=[-1,-1]; 
               stable = [fipie, fixe];
               stable2 = [stable2; stable];
            else
               detere1=[fipie,fixe];   
               detere=[detere; detere1];
            end



        end
    end

    Me1=sortrows(detere,2);



    display('INDETERMINACY REGION FOR FB-DITR')
    COMBe =Me1(obs1,:);  %Length of indeterminacy
    disp(COMBe)   
    save('ditre.mat', 'stable2', 'detere', 'COMBe');
end

%     figure;
%     hold on
% 
%         subplot(1,1,1)
%         plotmatrix(detere(:,1),detere(:,2))
%         title('Open economy with incomplete markets and FB-DITR');
%         xlabel('\phi_{\pi}');
%         ylabel('\phi_{x}');
% 
%     hold off
    

%==========================================================================
%     GENERATE PATCH SURFACES
%==========================================================================

colortheme1 = [.69  .69  .69];
colortheme2 = [.078  .745  .945];

xmax = 4;
ymax = 4;

% Patch diagram for DITR case
figure
    title('Open economy with incomplete markets and DITR');
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
    
    % Save this figure as .EPS and .FIG
    print('-depsc', 'oe-im-ditr')
    hgsave('oe-im-ditr')

% Patch diagram for FB-DITR case
figure
    title('Open economy with incomplete markets and FB-DITR');
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
    print('-depsc', 'oe-im-fbditr')
    hgsave('oe-im-fbditr')
    
    