%CALCBEDSHEARSTRESS.M
% This code takes the processed velocity information from the Vectrino
% Profiler or Vectrino and computes profiles of Reynolds, turbulent, and
% laminar stresses, from which bed shear stress can be computed. Before
% running this code, double click on the file name of the processed
% velocity file you would like to run to load it into the memory. 

% clear
clc
close all

%USER INPUT SECTION
nbins = sum(badbeams'<5&cell2mat(pct_rm)<0.5);%18;
bin_size = .002;
dist_above_bed = 6.4;
bin_range = 4.1:0.20:7.5; %4.05:0.2:6.45;
bin_range=bin_range(badbeams'<5&cell2mat(pct_rm)<0.5);
z = dist_above_bed-bin_range;
zmid = z(1:end-1)+0.1;
vpvp_mean=vpvp_mean(:,:,badbeams'<5&cell2mat(pct_rm)<0.5);
v = v(badbeams'<5&cell2mat(pct_rm)<0.5);

rho = 997; %kg/m^3
dyn_visc = 8.9e-4; %Pa*s
goodbins = intersect(find(badbeams<5), find(cell2mat(pct_rm)<0.5));

% figure(1), clf
% n = 1;
% labels = 'uvww';
% for ii = 1:4
%     for jj = 1:4
%         subplot(4,4,n)
%         plot(squeeze(vpvp_mean(ii,jj,:)), 1:13, 'ko')
%         title(sprintf('%s%s', labels(ii), labels(jj)))
%         n = n + 1;
%     end
% end

vmeans = NaN(nbins, 4);
uvmean = NaN(nbins,1);
uvpmean = NaN(1, nbins);
for ii = 1:length(goodbins)
    vmeans(ii, :) = mean(v{ii});
    latspeed = sqrt((v{ii}(:,1)).^2+(v{ii}(:,2)).^2);
    vertspeed = mean(v{ii}(:,3:4),2);
    uvpmean(ii) = mean((latspeed-mean(latspeed)).*(vertspeed-mean(vertspeed)));
end

speed = sqrt((vmeans(:,1)).^2+(vmeans(:,2)).^2+(mean([vmeans(:,3), vmeans(:,4)],2)).^2); %Better to pull this out.
       
% figure(2)
% for ii = 1:4
%     subplot(1,4,ii)
%     plot(-vmeans(:,ii), 1:13, 'ko')
% end

figure(3), clf

subplot(2,4,8)
plot(speed, z, 'ko')
title('Speed, m/s')

tau_TKE = 0.19*0.5*rho*(squeeze(vpvp_mean(1,1,:)+vpvp_mean(2,2,:))+transpose(mean(squeeze([vpvp_mean(3,3,:), vpvp_mean(4,4,:)]))));
figure(3), subplot(2,4,1), plot(tau_TKE, z, 'ko')
title('TKE')

tau_TKE_w = 0.9.*rho.*mean(squeeze([vpvp_mean(3,3,:), vpvp_mean(4,4,:)]));
figure(3), subplot(2,4,2), plot(tau_TKE_w, z, 'ko')
title('TKE w')

tau_Reynolds = -rho*uvpmean;
figure(3), subplot(2,4,3), plot(tau_Reynolds, z, 'ko')
title('Reynolds')

tau_laminar = dyn_visc*-diff(speed)./transpose((diff(z)/100));
figure(3), subplot(2,4,4), plot(tau_laminar, zmid, 'ko')
title('Laminar')

total_TKE = tau_laminar+tau_TKE(1:end-1)+diff(tau_TKE)/2;
subplot(2,4,5), plot(total_TKE, zmid, 'ko')

total_TKE_w = tau_laminar'+tau_TKE_w(1:end-1)+diff(tau_TKE_w)/2;
subplot(2,4,6), plot(total_TKE_w, zmid, 'ko')

total_Reynolds = tau_laminar'+tau_Reynolds(1:end-1)+diff(tau_Reynolds)/2;
subplot(2,4,7), plot(total_Reynolds, zmid, 'ko')

total_TKE(end-1)
total_TKE_w(end-3)
total_Reynolds(end-3)
speed(1)

% tau_TKE(end-5)
% tau_TKE_w(end-5)
% tau_Reynolds(end-5)