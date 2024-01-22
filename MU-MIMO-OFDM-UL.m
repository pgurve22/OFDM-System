clear all;
nSubC = 32;  
nCP = round(nSubC/4); 
niterations = 100;
nUEs = 4;
Ns=1;
Nr = 64;
No=1;
var_db = [0,-9.7,-19.2,-22.8,-30]; %Ped-A
%var_db = [0,-1,-9,-10,-15,-20]; %Veh-A
var = 10.^(var_db/10);
L=length(var_db);
SNR_db = -10:1:30;
SNR = 10.^(SNR_db/10);
beta_values = [0.749, 0.045, 0.246, 0.121, 0.125, 0.142, 0.635, 0.256, 0.673, 0.982, 0.1];

for combiner=1:3

for it=1:niterations
    
    h = zeros(Nr,nUEs,nSubC);
    for l=1:length(var_db)
        h(:,:,l)=sqrt(var(l)/2).*(randn(Nr,nUEs)+1j.*randn(Nr,nUEs));
    end

    for k=1:nUEs
        h(:,k,:) = sqrt(beta_values(k)).*(h(:,k,:));
    end

    for k=1:nUEs
        for m=1:nSubC
            Hm(:,k,m)=sum(reshape(h(:,k,:),Nr,nSubC).*exp(-1j*2*pi*(m-1).*((1:nSubC)-1)/(nSubC)),2);        
        end
    end

    for jj=1:length(SNR_db)
        [combiner,it,jj]
        if combiner==1 
            for m=1:nSubC
                W1(:,:,m)=Hm(:,:,m)';
            end
        elseif combiner==2
            for m=1:nSubC
                ch_zf=Hm(:,:,m).';
                W1(:,:,m)=pinv(ch_zf)';
            end
        elseif combiner==3
           for m=1:nSubC
                ch_mmse=Hm(:,:,m).';
                W1(:,:,m)=conj(ch_mmse)/(ch_mmse'*conj(ch_mmse)+((nUEs*No)/SNR(jj))*eye(Nr));
                %P1(:,:,m)=Hm(:,:,m)'/(Hm(:,:,m)*Hm(:,:,m)'+(nUEs*No/SNR(jj))*eye(nUEs));
            end
        end
    
        for k=1:nUEs
            for m=1:nSubC
                W(k,:,m) = W1(k,:,m)./(norm(W1(k,:,m))); 
            end
        end

        for m_dash=1:nSubC
            for k_dash=1:nUEs
                %Signal Power>>>>>>>>>>>>>>>>
                sgnpwr=abs(W(k_dash,:,m_dash)*Hm(:,k_dash,m))^2;
                num_term=sgnpwr;

                %Interference Power>>>>>>>>>>>>>>
                MUI=zeros(1,nUEs);
                for k=1:nUEs
                    if(k~=k_dash)
                        MUI(k)=abs(conj(W(k_dash,:,m_dash))*Hm(:,k,m_dash))^2;
                    end
                end
                dem_term=sum(MUI);

                SINR = (SNR(jj)*num_term)/(dem_term*SNR(jj)+1);
                SUM_rate(k_dash) = log2(1+SINR);
            end
            SUM_rate_avr(m_dash)=sum(SUM_rate);
        end
        SE(it,jj) = sum(SUM_rate_avr)/nSubC;
    end
end

SE_avr(combiner,:) = sum(SE,1)/niterations;
end

plot(SNR_db,SE_avr(1,:),'LineWidth', 1.8, 'Marker', 'square', 'MarkerSize', 8 , 'MarkerFaceColor','g');
hold on;
plot(SNR_db,SE_avr(2,:),'LineWidth', 1.8, 'Marker', 'square', 'MarkerSize', 8 , 'MarkerFaceColor','b');
plot(SNR_db,SE_avr(3,:),'LineWidth', 1.8, 'Marker', 'square', 'MarkerSize', 8 , 'MarkerFaceColor','c');
legend('MRC','ZF','MMSE','Location','northwest');
xlabel("SNR");
ylabel("Uplink Sum Rate");
grid on;






