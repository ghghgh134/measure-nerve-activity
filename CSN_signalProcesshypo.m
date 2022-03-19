clear
load SHR01lowo2.mat
CSNA1 = data;
amp = 50000;
CSNA1 = CSNA1/amp*10^6;% baseline under 100% oxygen
T1 = (0:0.1:0.1*(length(CSNA1)-1));% baseline under 100% oxygen
T1 = T1*0.005;
T1 = T1(1:120001);
CSNA1 = CSNA1(1:120001);

fs = 2000;
CSNA1 = highpass(CSNA1,50,fs);

load SHR01PElowo2.mat
CSNA2 = data;
amp = 50000;
CSNA2 = CSNA2/amp*10^6;% PE infuse 3min 100% oxygen
T2 = (0:0.1:0.1*(length(CSNA2)-1));
T2 = T2*0.005;
T2 = T2(1:300001);
CSNA2 = CSNA2(1:300001);

fs = 2000;
CSNA2 = highpass(CSNA2,50,fs);

load SHR01SNPlowo2.mat
CSNA3 = data;
amp = 50000;
CSNA3 = CSNA3/amp*10^6;% SNP infuse 2.5min 100% oxygen
T3 = (0:0.1:0.1*(length(CSNA3)-1));
T3 = T3*0.005;
T3 = T3(1:300001);
CSNA3 = CSNA3(1:300001);

fs = 2000;
CSNA3 = highpass(CSNA3,50,fs);

CSNA = {CSNA1 CSNA2 CSNA3};
To2 = {T1 T2 T3};
figure(1);clf;hold on;
subplot(2,3,4)
plot(T1,CSNA1);
ylim([-100 100]);

subplot(2,3,5)
plot(T2,CSNA2);
ylim([-100 100]);

subplot(2,3,6)
plot(T3,CSNA3);
ylim([-100 100]);
%%
load SHR01.txt
Tbp = SHR01(:,1);
BP = SHR01(:,2);
% Creating smoothing window for smooth pressure
dt = 0.02;
tw = -2:dt:2;
W = exp(-(tw.^2)/0.001); W = W./sum(W);
SmoothPressureextract = conv(BP,W);  BP = SmoothPressureextract(51:length(BP)+50);
Tmarker = xlsread('SHR01.xlsx',1);
Start = Tmarker(:,1); 
End = Tmarker(:,2); 
BPseq = cell(length(Start),1);
Tseq = cell(length(Start),1);
for i = 1:length(Start)
    p = find((Tbp <= End(i))&(Tbp >= Start(i)));
    BPseq{i} = BP(p);
    Tseq{i} = (0:0.002:(length(BPseq{i})-1)*0.002)';
end
BPseq = BPseq(4:6);
Tseq = Tseq(4:6);
figure(1);
subplot(2,3,1);
plot(Tseq{1},BPseq{1});hold on;
subplot(2,3,2);
plot(Tseq{2},BPseq{2});hold on;
subplot(2,3,3);
plot(Tseq{3},BPseq{3});hold on;
%%
Tnew = cell(1,3);
DP = cell(1,3);
SP = cell(1,3);
MAP = cell(1,3);
PP = cell(1,3);
HR = cell(1,3);
RR = cell(1,3);
deltaS = cell(1,3);
INA = cell(1,3);
M = cell(1,3);
P = cell(1,3);
ina = cell(1,3);

for i = 1 : length(Tseq)
    Tnew{i} = (Tseq{i}(1):0.0001:Tseq{i}(end));
    % Diastolic BP
    [dp,dp_locs] = findpeaks(-BPseq{i},'MinPeakProminence',20); 
    d = griddedInterpolant(Tseq{i}(dp_locs),-dp,'pchip'); 
    DP{i} = d(Tnew{i}); 
    
    % Systolic BP
    [sp,sp_locs] = findpeaks(BPseq{i},'MinPeakProminence',20); 
    s = griddedInterpolant(Tseq{i}(sp_locs),sp,'pchip'); 
    SP{i} = s(Tnew{i});
    
    % Mean arterial BP and Pulse BP
    for j = 1:length(dp_locs)-1
     int = dp_locs(j):dp_locs(j+1);
     M{i}(j) = trapz(Tseq{i}(int),BPseq{i}(int))/(Tseq{i}(int(end)) - Tseq{i}(int(1))); 
     P{i}(j) = max(BPseq{i}(int))-min(BPseq{i}(int));

    end
    x = find(~isnan(M{i}));
    t = Tseq{i}(dp_locs(x)); 
    m = griddedInterpolant(t,M{i}(x),'pchip');
    p = griddedInterpolant(t,P{i}(x),'pchip');
    MAP{i} = m(Tnew{i}); 
    PP{i} = p(Tnew{i});
    % Heart rate
    HR_raw = 1./diff(Tseq{i}(dp_locs))*60; 
    HR_raw = [HR_raw(1); HR_raw]; 
    hsmooth = smooth(HR_raw,0.005,'rloess');
    h = griddedInterpolant(Tseq{i}(dp_locs),hsmooth,'pchip'); 
    HR{i} = h(Tnew{i}); 
    
    % RR interval 
    RR{i} = 1./HR{i}*60; 


    [ns,ns_locs] = findpeaks(-CSNA{i},'MinPeakProminence',2); 
    d = griddedInterpolant(To2{i}(ns_locs),-ns,'pchip'); 
    NS = d(Tnew{i}); 
    
    [ps,ps_locs] = findpeaks(CSNA{i},'MinPeakProminence',2); 
    s = griddedInterpolant(To2{i}(ps_locs),ps,'pchip'); 
    PS = s(Tnew{i});
    deltaS{i} = PS-NS;
    % Mean sigal

    for j = 1:length(dp_locs)-1
     int = find(Tnew{i} > Tseq{i}(dp_locs(j)) & Tnew{i} < Tseq{i}(dp_locs(j+1)));
     ina{i}(j) = trapz(Tnew{i}(int),deltaS{i}(int))/(Tnew{i}(int(end)) - Tnew{i}(int(1))); 
    end
    x = find(~isnan(ina{i}));
    t = Tseq{i}(dp_locs(x)); 
    in = griddedInterpolant(t,ina{i}(x),'pchip');
    INA{i} =in(Tnew{i}); 
end
figure(2);clf; hold on;
subplot(3,3,1);
plot(Tnew{1},PP{1});
subplot(3,3,4);
plot(Tnew{1},MAP{1});
subplot(3,3,7);
plot(Tnew{1},INA{1});
subplot(3,3,2);
plot(Tnew{2},PP{2});
subplot(3,3,5);
plot(Tnew{2},MAP{2});
subplot(3,3,8);
plot(Tnew{2},INA{2});
subplot(3,3,3);
plot(Tnew{3},PP{3});
subplot(3,3,6);
plot(Tnew{3},MAP{3});
subplot(3,3,9);
plot(Tnew{3},INA{3});
figure(3);clf;
plot(T1,CSNA1);hold on;
plot(Tnew{1},deltaS{1});hold on;
figure(4);clf;
plot(T2,CSNA2);hold on;
plot(Tnew{2},deltaS{2});hold on;
figure(5);clf;
plot(T3,CSNA3);hold on;
plot(Tnew{3},deltaS{3});hold on;
M = [M{1}' ; M{2}' ; M{3}' ];
P = [P{1}' ; P{2}' ; P{3}' ];
ina = [ina{1}'; ina{2}' ;ina{3}'];

%%


nM = [min(M):1:max(M)];
nina = interp1(M,ina,nM,'pchip');
D1 = [nM' nina'];
D1 = sortrows(D1);


% Fit model to data.
figure(6);clf
[~,fitresult1] = sigm_fit(D1(:,1), D1(:,2) );
scatter(nM,nina,'filled','Markeredgecolor',[0.7 0.7 0.7],'sizedata',3,'Markerfacecolor',[0.7 0.7 0.7]);hold on;
plot(D1(:,1),fitresult1.ypred,'linewidth',1.5);hold off;
save lowo2.mat nM fitresult1 nina

nP = [min(P):1:max(P)];
nina = interp1(P,ina,nP,'pchip');
D2 = [nP' nina'];
D2 = sortrows(D2);

% Fit model to data.
figure(7);clf
[~,fitresult2] = sigm_fit(D2(:,1), D2(:,2) );
scatter(nP,nina,'filled','Markeredgecolor',[0.7 0.7 0.7],'sizedata',3,'Markerfacecolor',[0.7 0.7 0.7]);hold on;
plot(D2(:,1),fitresult2.ypred,'linewidth',1.5);hold off;