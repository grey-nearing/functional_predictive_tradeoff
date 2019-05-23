%clear all; close all; clc;
restoredefaultpath; addpath('matlab_tools');

%% *** Dimensions ********************************************************

if 1

% data directory
datadir = './data/raw_data/';

% number of info bins
Nbins = 30;

% number of files
Nfiles = 118;

% training fraction
trnfrac = 0.7;
Nepochs = 2000;

% lags
lags = [1];%,2,4,6];%,8,12];
Nlags = length(lags);

% input/output/obs indexes
Iu = [7,11,12,16,17,31,32]; Du = length(Iu);
Iy = [1,4,8,13];            Dz = length(Iy);
Iz = [2,5,9,14];
%Ix = [1,4,7,8,11,12,13,16,17,18:30,31,32]x];

%% *** Gather Data *******************************************************
    
% init storage
TE = zeros(Du+Dz,Dz,Nlags,Nfiles,3)/0;
MI = zeros(Du+Dz,Dz,Nlags,Nfiles,3)/0;
Hx = zeros(Du+Dz,Dz,Nlags,Nfiles,3)/0;
Hy = zeros(Du+Dz,Dz,Nlags,Nfiles,3)/0;

% loop thgouth models
for f = 1:Nfiles
    
    % load raw data
    fname = strcat(datadir,num2str(f),'.txt');
    data = load(fname);
    data(data<=-9990) = 0/0;
    Nt = size(data,1);
    
    % extract raw data
    U = data(:,Iu);
    Y = data(:,Iy);
    Z = data(:,Iz);
    
    % check that the obs are everywhere constant
    if f == 1; Ukeep = U; Zkeep = Z; end
    I = find(~isnan(U(:)));     if any(U(I)~=Ukeep(I)); continue; end
    I = find(~isnan(Z(:)));     if any(Z(I)~=Zkeep(I)); continue; end
    I = find(~isnan(Ukeep(:))); if any(U(I)~=Ukeep(I)); continue; end
    I = find(~isnan(Zkeep(:))); if any(Z(I)~=Zkeep(I)); continue; end
    
    % pull good data
    XXmod = [U,Y];
    XXobs = [U,Z];
    YYmod = [Y];
    YYobs = [Z];
    
    % dimensions
    Dx = Du+Dz; assert(size(XXobs,2)==Dx); assert(size(XXmod,2)==Dx);
    Dy = Dz;    assert(size(YYobs,2)==Dy); assert(size(YYmod,2)==Dy);
    
    % loop through lags
    for l = 1:Nlags
        
        % screen report
        tic;
        fprintf('model file %d/%d -- lag %d/%d ...',f,Nfiles,l,Nlags);
        
        % lag data
        clear Xo Xm
        for x = 1:Dx
            Xo(:,x) = window_average(XXobs(:,x),lags(l));
            Xm(:,x) = window_average(XXmod(:,x),lags(l));
        end
        
        clear Yo Ym
        for y = 1:Dy
            Yo(:,y) = window_average(YYobs(:,y),lags(l));
            Ym(:,y) = window_average(YYmod(:,y),lags(l));
        end
        
        % check dimensions
        N = size(Yo,1);
        assert(N==size(Ym,1));
        assert(N==size(Xo,1));
        assert(N==size(Xm,1));
        
        % lag data
        XoL = Xo(2:end,:);
        XmL = Xm(2:end,:);
        YoL = Yo(1:end-1,:);
        YmL = Ym(1:end-1,:);
        
        % deal with grandma
        I = find(any(isnan(XoL'))); XoL(I,:) = []; YoL(I,:) = []; XmL(I,:) = []; YmL(I,:) = [];
        I = find(any(isnan(XmL'))); XoL(I,:) = []; YoL(I,:) = []; XmL(I,:) = []; YmL(I,:) = [];
        
        % make sure we have dealt with missing values
        assert(isempty(find(isnan(XoL),1,'first')));
        assert(isempty(find(isnan(XmL),1,'first')));
        
        for y = 1:Dy
            
            % pull data for this pathway
            XToL = XoL;
            XTmL = XmL;
            YSoL = YoL(:,y);
            YSmL = YmL(:,y);
            
            % deal with pathway-specific grandmas
            I = find(isnan(YSoL)); YSoL(I) = []; XToL(I,:) = []; YSmL(I) = []; XTmL(I,:) = [];
            I = find(isnan(YSmL)); YSoL(I) = []; XToL(I,:) = []; YSmL(I) = []; XTmL(I,:) = [];
            Ntot = length(YSoL);
            
            % make sure we have dealt with missing values
            assert(isempty(find(isnan(XToL),1,'first')));
            assert(isempty(find(isnan(XTmL),1,'first')));
            assert(isempty(find(isnan(YSoL),1,'first')));
            assert(isempty(find(isnan(YSmL),1,'first')));
            
            % regression models
            Ntrn = round(trnfrac*Ntot);
            if f == 1;
                Itrn = randperm(Ntot,Ntrn);
                Itst = setdiff(1:Ntot,Itrn);
                [regModel(l,y),~,RSoL] = ann_experiment(XToL,YSoL,Itrn,Itst,'trainscg',Nepochs,0);
                Rkeep{l,y} = RSoL;
            else
                [~,RSoL] = pred_ann(XToL,YSoL,regModel(l,y),1:Ntrn,[]);
                same = find(any(Rkeep{l,y}~=RSoL));
                assert(isempty(same));
            end
            
            % info bins
            Bmin = min([min(YSmL),min(YSoL),min(RSoL)])-1e-6;
            Bmax = max([max(YSmL),max(YSoL),max(RSoL)])+1e-6;
            db =(max(YoL(:,y))-min(YoL(:,y)))/Nbins;
            By = Bmin:db:(Bmax+db);
            
            % loop through pathways
            for x = 1:Du+Dy
                
                % special case when targeting self
                YToL = XToL(:,y+Du);
                YTmL = XTmL(:,y+Du);
                if x == Du+y
                    YToL = rand(size(YToL)) * (max(YToL)-min(YToL)) + min(YToL);
                    YTmL = rand(size(YTmL)) * (max(YTmL)-min(YTmL)) + min(YTmL);
                end
                
                % info bins
                Bmin = min([min(XmL(:,x)),min(XoL(:,x))])-1e-6;
                Bmax = max([max(XmL(:,x)),max(XoL(:,x))])+1e-6;
                db =(max(XoL(:,x))-min(XoL(:,x)))/Nbins;
                Bx = Bmin:db:(Bmax+db);
                
                % transfer entropies
                [TE(x,y,l,f,1),Hx(x,y,l,f,1),Hy(x,y,l,f,1)] = transfer_entropy(YSmL,XTmL(:,x),YTmL,Bx,By); % model pathway
                [TE(x,y,l,f,2),Hx(x,y,l,f,2),Hy(x,y,l,f,2)] = transfer_entropy(YSoL,XToL(:,x),YToL,Bx,By); % observation pathway
                [TE(x,y,l,f,3),Hx(x,y,l,f,3),Hy(x,y,l,f,3)] = transfer_entropy(RSoL,XToL(:,x),YToL,Bx,By); % regression pathway
                
                % mutual informations
                [MI(x,y,l,f,1),~,~] = mutual_info(XTmL(:,x),YSmL,Bx,By); % model pathway
                [MI(x,y,l,f,2),~,~] = mutual_info(XToL(:,x),YSoL,Bx,By); % observation pathway
                [MI(x,y,l,f,3),~,~] = mutual_info(XToL(:,x),RSoL,Bx,By); % regression pathway
                
                % calculate statistics
                [iyz,hz,~] = mutual_info(YSoL,YSmL,By,By);
                mmi(f,y) = iyz/hz;
                [iyz,hz,~] = mutual_info(RSoL,YSmL,By,By);
                rmi(f,y) = iyz/hz;
                mae(f,y) = mean(abs((YSoL-YSmL)./YSoL));
                rae(f,y) = mean(abs((RSoL-YSmL)./RSoL));
                mse(f,y) = sqrt(mean((YSoL-YSmL).^2));
                rse(f,y) = sqrt(mean((RSoL-YSmL).^2));
                
            end % y
        end % x
        
        % screen report
        t = toc;
        fprintf('finished = %f \n',t);
        
    end % lags
end % model files

%% --- Save Results -------------------------------------------------------

% save results
fname = 'dpn_ann_analysis_results.mat';
save(fname);

else
    fname = 'dpn_ann_analysis_results.mat';
    load(fname);
end

%% --- Plot Results -------------------------------------------------------

close all
fignum = 0;

% Variable Names
Unames = [{'LAI'},{'P'},{'Rg'},{'T'},{'U'},{'VPD'},{'ZEN'}];
Ynames = [{'FCO_2'},{'H'},{'LE'},{'Rn'}];
Pnames = [{'Cd'},{'CO2'},{'mslope'},{'P'},{'Ta'},{'Vcmax25'}];
Np = [24,20,24,11,15,24];
Anames = [Unames,Ynames];

% markers
markers = ['o','s','d','h','p','^','v','<','>','o','s'];

% colors
clear colors;
figure(1); h = plot(rand(7));
for i = 1:7
    colors(i,:) = h(i).Color;
end; close(1); clear h;
r = linspace(0,1,Du+Dz);
for i = 8:Du+Dz
    p = randperm(11,3);
    colors(i,:) = [r(p(1)),r(p(2)),r(p(3))];
end

fignum = 1;
for p = [1,3,6]

    fignum = fignum+1;
    figure(fignum); close(fignum); figure(fignum);
    set(gcf,'color','w','position',[1072           1         329         804]);
    splot = 0;
    
    for y = 1:Dz-1
        splot = splot+1;
        subplot(Dz-1,1,splot);
        
        idex = sum(Np(1:p-1))+1:sum(Np(1:p));
        
        for x = 1:Du
            TEdiffO = squeeze((TE(x,y,1,:,1) - TE(x,y,1,:,2))./TE(x,y,1,:,1));
            TEdiffR = squeeze((TE(x,y,1,:,1) - TE(x,y,1,:,3))./TE(x,y,1,:,1));
            
            plot(TEdiffO(idex),1-mmi(idex,y),...
                'linestyle','none',...
                'color','k',...
                'marker',markers(x),...
                'markersize',10,...
                'markerfacecolor',colors(x,:));
            hold on;
        end
        
        set(gca,'fontsize',12)
        xl = xlim; yl = ylim; %yl(1) = 0;
        plot([0,0],[yl(1),yl(2)],'--k','linewidth',2);
        grid on
        if splot == 3
            leg = legend(Unames,'location','best'); 
%             set(leg,'position',[0.80547      0.80037      0.18997      0.12189]);
        end
        axis([xl,yl]);
        
%         xlabel('\Delta Transfer Entropy (%)','fontsize',18);
%         ylabel('Norm. Missing MI [nats/nats]','fontsize',18);
        xlabel(strcat({'Functional Performance A_f (y_i -> '},Ynames{y},')'),'fontsize',14);
        ylabel(strcat('Predictive Performance A_p (',Ynames{y},')'),'fontsize',14);
        titstr = strcat({'Varying '''},Pnames{p},{''' to Estimate '},Ynames{y});
        title(titstr,'fontsize',16);
        
    end

    % save
    fname = strcat('./figures/Fig3 - ParetoTradeoffs_',Pnames{p},'.png');
    figure(fignum);
    set(gcf,'PaperPositionMode','auto')
    saveas(gcf,fname);
    
end



%% --- Close-up Plot -------------------------------------------------------

% initialize figure
fignum = fignum+1;
figure(fignum); close(fignum); figure(fignum);
set(gcf,'color','w','position',[1120         362         824         383]);
splot = 0;

for y = 1
    for p = 3
        splot = splot+1;
        subplot(1,1,splot);
        
        idex = sum(Np(1:p-1))+1:sum(Np(1:p));
        
        for x = 1:Du
            TEdiffO = squeeze((TE(x,y,1,:,1) - TE(x,y,1,:,2)));%./TE(x,y,1,:,1));
            TEdiffR = squeeze((TE(x,y,1,:,1) - TE(x,y,1,:,3)));%./TE(x,y,1,:,1));
            
            plot(TEdiffO(idex),1-mmi(idex,y),...
                'linestyle','none',...
                'color','k',...
                'marker',markers(x),...
                'markersize',10,...
                'markerfacecolor',colors(x,:));
            hold on;
        end
        
        set(gca,'fontsize',18)
        xl = xlim; yl = ylim;
        plot([0,0],[yl(1),yl(2)],'--k','linewidth',2);
        grid on
%         set(gca,'xticklabels',[],'yticklabels',[]);
%         if y == 1
%             axis([-0.8,0.2,1-0.399,1-0.392]);
%         elseif y == 2
%             axis([0,70,0.24,0.265]);
%         elseif y == 3
%             axis([-150,50,0.44,0.455]);
%         end
        
%        xlabel('\Delta Transfer Entropy (%)','fontsize',18);
%        ylabel('Normalized MI [nats/nats]','fontsize',18);
%        titstr = strcat({'Varying '''},Pnames{p},{''' to Estimate '},Ynames{y});
%        title(titstr,'fontsize',22);
        
    end
end

% save
figure(fignum);
set(gcf,'PaperPositionMode','auto')
saveas(gcf,'./figures/Fig3 - ParetoTradeoffs(subplot).png');

%% --- Figure 4 -----------------------------------------------------------

% initialize figure
fignum = fignum+1;
figure(fignum); close(fignum); figure(fignum);
set(gcf,'color','w','position',[793   255   329   533]);
splot = 0;

p = 3;
idex = sum(Np(1:p-1))+1:sum(Np(1:p));
y = 1;

for x = 1:Du
    TEdiffO = squeeze((TE(x,y,1,:,1) - TE(x,y,1,:,2))./TE(x,y,1,:,1));
    TEdiffR = squeeze((TE(x,y,1,:,1) - TE(x,y,1,:,3))./TE(x,y,1,:,1));
    
    subplot(2,1,1);
    plot(TEdiffO(idex),1-mmi(idex,y),...
        'linestyle','none',...
        'color','k',...
        'marker',markers(x),...
        'markersize',10,...
        'markerfacecolor',colors(x,:));
    hold on;
    
    subplot(2,1,2);
    plot(TEdiffR(idex),1-mmi(idex,y),...
        'linestyle','none',...
        'color','k',...
        'marker',markers(x),...
        'markersize',10,...
        'markerfacecolor',colors(x,:));
    hold on;
end

for s = 1:2
    subplot(2,1,s);
    set(gca,'fontsize',12)
    xl = xlim; yl = ylim;
    plot([0,0],[yl(1),yl(2)],'--k','linewidth',2);
    grid on
%     axis([-0.05,0.03,1-0.399,1-0.392]);
%     axis([-0.08,0.08,0.54,0.66]);
    axis([-0.8,0.2,1-0.399,1-0.392]);
    xlabel(strcat({'Functional Performance A_f (y_i -> '},Ynames{y},')'),'fontsize',14);
    ylabel(strcat('Predictive Performance A_p (',Ynames{y},')'),'fontsize',14);
end

subplot(2,1,1);
titstr = 'Raw Obs: A_f = [I_m_o_d - I_o_b_s]/I_m_o_d';
title(titstr,'fontsize',16);

subplot(2,1,2);
titstr = 'ANN: A_f = [I_m_o_d-I_r_e_g]/I_m_o_d';
title(titstr,'fontsize',16);
legend(Unames,'location','sw');

% save
figure(fignum);
set(gcf,'PaperPositionMode','auto')
saveas(gcf,'./figures/FigA1 - obsVann.png');

%% --- END PROGRAM --------------------------------------------------------







