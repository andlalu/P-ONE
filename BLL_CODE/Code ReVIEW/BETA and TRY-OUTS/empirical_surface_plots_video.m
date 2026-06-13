clc; clear; 
load('prelucrat_data9610v2.mat');
%%%%%%%%%%%%
% Gaussian kernel smoothing to prepare standardized surface plots:
        newdata_graph=[];
        %Set-up the implied volatility surface gridpoints for smoothing and plotting:
        std_T_graph=(0.1:0.05:1);
        std_K_graph=(0.8:0.01:1.2);
        gridpoints_graph = [reshape(repmat(std_T_graph,[length(std_K_graph),1]),[length(std_K_graph)*length(std_T_graph),1]),reshape(repmat(std_K_graph',[1,length(std_T_graph)]),[length(std_K_graph)*length(std_T_graph),1])];
        %Apply the kernel smoothing procedure for each sample period 'i' 
%         for i=1:max(fulldata(:,1))'
for i=655:680
             index_i=(fulldata(:,1)==i);
             if sum(index_i)>0 %check if there are prices quoted on date 'i' 
             forward_sel=unique(fulldata(fulldata(:,1)==i,[2,6]),'rows'); %forward rate curve
             forwards_graph=spline(forward_sel(:,1),forward_sel(:,2),std_T_graph); %interpolate forward rate curve  
             forwards_graph=repmat(forwards_graph,[length(std_K_graph),1]);
             forwards_graph=forwards_graph(:);
             newdata_graph=[newdata_graph; [i.*ones([length(gridpoints_graph),1]),gridpoints_graph,kernel_vol(gridpoints_graph, [fulldata(index_i,2),fulldata(index_i,3)./fulldata(index_i,6),fulldata(index_i,7)])]];
             end
        end
%%%%%%%%%%%%
% Implied volatility surface plots: 
    
x=663;
    for t=x+4
    hold on;
    volplot(std_T_graph,std_K_graph,reshape(newdata_graph(newdata_graph(:,1)==t,4),[length(std_K_graph),length(std_T_graph)]))
    hold off;
    end

% dynamics(std_T_graph,std_K_graph,...
%     reshape(newdata_graph(newdata_graph(:,1)==x,4),[length(std_K_graph),length(std_T_graph)]),...
%     reshape(newdata_graph(newdata_graph(:,1)==x+8,4),[length(std_K_graph),length(std_T_graph)])),...
%     %reshape(newdata_graph(newdata_graph(:,1)==x,4),[length(std_K_graph),length(std_T_graph)])),...
    %reshape(newdata_graph(newdata_graph(:,1)==x+12,4),[length(std_K_graph),length(std_T_graph)]),...
    %reshape(newdata_graph(newdata_graph(:,1)==x+16,4),[length(std_K_graph),length(std_T_graph)]));


%     surf(reshape(newdata_graph(newdata_graph(:,1)==t,4),[length(std_K_graph),length(std_T_graph)]))
    %surf(std_T_graph',std_K_graph',reshape(newdata_graph(newdata_graph(:,1)==t,4),[length(std_T_graph),length(std_K_graph)]))%,'EdgeColor',((25-t)/25).*[0.85 0.85 0.85],'FaceColor',((25-t)/25).*[0.9 0.9 0.9]);
% %% Set up the movie.
% writerObj = VideoWriter('out.avi'); % Name it.
% writerObj.FrameRate = 60; % How many frames per second.
% open(writerObj); 
%  
% for i=1:size(y)      
%     % We just use pause but pretend you have some really complicated thing here...
%     pause(0.1);
%     figure(fId); % Makes sure you use your desired frame.
%     plot(x(i),y(i),'or');
%  
%     %if mod(i,4)==0, % Uncomment to take 1 out of every 4 frames.
%         frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
%         writeVideo(writerObj, frame);
%     %end
%  
% end
% hold off
% close(writerObj); % Saves the movie.