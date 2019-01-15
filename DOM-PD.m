WLset=284:460;
nWL=177;
lambda=270:700;
[~,WLindex,~]=intersect(lambda,WLset); % indices for astar values
%
% .  initial values
%
CDOC=zeros(700,3);

pdCDOC=zeros(3,4);
% AQYvec=exp(-0.02767*(WLset-284));
AQYvec3=exp((-0.0210)*(WLset-284));
AQYvec2=exp((-0.0495)*(WLset-284));
AQYvec1=exp((-0.0328)*(WLset-284));
%
% . iterate time steps of one minute for 12 hours
ICM_input=importdata('/Users/bclark/Desktop/test_WC_PD/wqm_kei_photoDEG_v8.csv');
load('Incubation_data_v8.mat');
lambda=ICM_input(:,1);
astar_cdom1 = ICM_input(:,3);
astar_cdom2 = ICM_input(:,4);
astar_cdom3 = ICM_input(:,5);

astarmdl=[astar_cdom1 astar_cdom2 astar_cdom3];
mkdir('figures');

AQYnew_mol=[];
%%
%  [naqy,~]=size(good_APQY);

for iapqy =1 :1;%naqy;
    
%  APQY_in=good_APQY(iapqy,:)./12.*1000;
%  APQY_in=median(good_APQY,2)./12*1000/10;
    
%     
% % 
% %    % Marsh low tide only
    AQYnew_mol(3,1,:)=5.74.*AQYvec3;%APQY_in(2);
    
    AQYnew_mol(3,2,:)=0.0;%APQY_in(1);
    AQYnew_mol(3,3,:)=0.111*0.5.*AQYvec3;%APQY_in(3);
    AQYnew_mol(3,4,:)=0.111*0.5.*AQYvec3;%APQY_in(4);
    AQYnew_mol(2,1,:)=56.05.*AQYvec2;%APQY_in(5);
    AQYnew_mol(2,2,:)=0.0;
    AQYnew_mol(2,3,:)=3.23E-8*0.5.*AQYvec2;%APQY_in(6);
    AQYnew_mol(2,4,:)=3.23E-8*0.5.*AQYvec3;%APQY_in(7);
    AQYnew_mol(1,3,:)=1.665*0.5.*AQYvec1;%APQY_in(8);
    AQYnew_mol(1,4,:)=1.665*0.5.*AQYvec1;%APQY_in(9);
    % Marsh low tide only
    Iz0=Expirrtest./1000; % cut out the section of spectra we want and convert to W
 
    id_300 = find(lambda == 300); % find the 300 nm wavelentght
    s1=find(lambda >=275 & lambda <=295);   % the two ares to pull out our slopes from
    s2 = find(lambda >=350& lambda <=400);
%           figure;
          CDOC_out=[];
         
  %  uncomment here if want to use average marsh low tide or estuarine conditions     
%      CDOC_in=CDOC_avg;
 [ndom,~]=size(CDOC_in);
    for idom=1:ndom;
        
        CDOC(1,:)=CDOC_in(idom,:);  %Average Marsh CDOC from Clark et al.
        pdCDOC_out=[];
        for ts=1:700
            %  CALCULATE absorbed irradiance
            %
            absmdl=repmat(CDOC(ts,:),nWL,1).*astarmdl(WLindex,:);
            
            Izmdl_in=repmat(Iz0.*(1-exp(-sum(absmdl,2)*.04))./sum(absmdl,2),1,3).*absmdl;%
            %

            Izmdl=Izmdl_in.*WLset'.*5.03E15./6.02249E23;
      
            
            for i=1:3 % source pools
                
                for j=1:4 % product pools j=3: mineralization (0) . j=4:non-colored DOC
                    
                    pdCDOC(i,j)=sum(Izmdl(:,i).*squeeze(AQYnew_mol(i,j,:)))./.04;%12 e-3 converts mmol to g to give result in g/m3=mg/L

                end
            end
            
            Izmdl_out(:,ts,idom)=sum(Izmdl,2);
            %
            % calculate change in each pool rate is mg/L/min
            %
            pdCDOC_out=cat(3,pdCDOC_out,pdCDOC)    ;
            pdCDOC=pdCDOC.*12E-3*864; %convert from mmol C m^-3 s^-1 to g C m^-3 timestep^-1 (time step is 0.01 days)
            
          
            dCDOC(ts,1,idom)=pdCDOC(3,1)+pdCDOC(2,1)-pdCDOC(1,3)-pdCDOC(1,4);
            dCDOC(ts,2,idom)=pdCDOC(3,2)-pdCDOC(2,1)-pdCDOC(2,3)-pdCDOC(2,4);
            dCDOC(ts,3,idom)=-pdCDOC(3,1)-pdCDOC(3,2)-pdCDOC(3,3)-pdCDOC(3,4);
            %
            % adjust pool sizes for this time step
            %
            if(ts<701)
                CDOC(ts+1,:)=CDOC(ts,:)+dCDOC(ts,:,idom);
            end
            
        end
        pdCDOC_out2(:,:,:,idom)=pdCDOC_out;
        CDOC_out(:,:,idom)=CDOC;
    end

    
%     SR_obs=[2.16 1.81 2.12 1.96	1.78 1.99 1.45	1.72 1.55 1.89	2.15 1.97 1.59 1.72	1.64];
%     
%     S2_obs=[-0.013279837	-0.013583881	-0.012662487	-0.016360128	-0.01382052	-0.013781571	-0.014824009	-0.014351294	-0.014609292	-0.014522053	-0.0139961	-0.0150	-0.0145 -0.0148 -0.0153];
%     
%     S1_obs=[-0.0286	-0.0246	-0.0269	-0.0320	-0.0246	-0.0273	-0.0214	-0.0246	-0.0226	-0.0275	-0.0301...
%         -0.029	-0.0231	-0.0254	-0.0250];

% observed final absorbance at 300 nm
  a300_obs=[3.441,8.360,5.423,3.097,8.861,3.983,15.222,11.453,20.533,9.692,2.806,3.904,16.103,10.111,9.448];
    clear big_abs_total big_s1 big_s2 big_sr absorbance
    for i =1:ndom;
       
        a300_out(:,i)=CDOC_out(end,:,i).*astarmdl(31,:);
        
        %       disp(['calculating absorbance for model run ',num2str(j)])
        % calculated absorption in different wave bands
        for ilamb = 1 : length(s1)
            as1(:,1,ilamb) = CDOC_out(:,1,i).*astarmdl(s1(ilamb),1);
            as1(:,2,ilamb) = CDOC_out(:,2,i).*astarmdl(s1(ilamb),2);
            as1(:,3,ilamb) = CDOC_out(:,3,i).*astarmdl(s1(ilamb),3);
            %         disp(['calculating absorbance for wavelength ',num2str(ilamb)])
        end
        
        for ilamb = 1 : length(s2)
            
            as2(:,1,ilamb) = CDOC_out(:,1,i).*astarmdl(s2(ilamb),1);
            as2(:,2,ilamb) = CDOC_out(:,2,i).*astarmdl(s2(ilamb),2);
            as2(:,3,ilamb) = CDOC_out(:,3,i).*astarmdl(s2(ilamb),3);
            %         disp(['calculating absorbance for wavelength ',num2str(ilamb)])
        end
        % end
        
        total_s1  = sum(as1,2);
        total_s2  = sum(as2,2);
        
        % now find the SR and slope regions
        
        [times,ids]=size(total_s1);
        
        for t = 1 : times
            %         disp(['Find the slopes for our regions at time ', num2str(t)])
            
            a=regress(log(total_s1(t,:)'),[ones(length(s1),1) [275:295]']);   % REGRESSION TO GET THE SLOPE FROM PAT NEALES CARY PROCESSING SCRIPT
            big_s1(t,i)=a(2);
            
            a=regress(log(total_s2(t,:)'),[ones(length(s2),1) [350:400]']);   % REGRESSION TO GET THE SLOPE FROM PAT NEALES CARY PROCESSING SCRIPT
            big_s2(t,i)=a(2);
            
        end
        
       for id = 1 :3;
        a300_ts(:,i,id)=CDOC_out(:,id,i).*astarmdl(31,id);
       end 
    big_SR = big_s1./big_s2;
 
    end

figure

a300_total=sum(a300_out,1);
[MEF(iapqy),r(iapqy),WMS(iapqy),RMSE(iapqy),MPE(iapqy),RI(iapqy),axis_min,axis_max]=model_statistics(a300_obs,a300_total);

for i = 1:ndom
    
    text(a300_obs(i),a300_total(i),num2str(i),'color','k');
    
    hold on
    % axis([min([a300_obs' a300_total']) max([a300_obs' a300_total']) min([a300_obs' a300_total']) max([a300_obs' a300_total'])]);
end
axis([axis_min axis_max axis_min axis_max]);
xlabel('Obs a300');ylabel('Pred a300');

%%
saveas(gcf,['figures/Scatter_out_' num2str(iapqy) '.png']);
saveas(gcf,['figures/Scatter_out_' num2str(iapqy) '.fig']);
saveas(gcf,['figures/Scatter_out_' num2str(iapqy) '.eps']);

end





























