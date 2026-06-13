%%%% 01.05.2015 V2.0 Final: for any AJD model.
function [ deriv_data, states_data, TT_marks, coeff_mat ] = sim_deriv_data(X0,N,dt,type,theta,varargin)
%sim_deriv_data Simulates a data-set of call derivative prices


%% Simulate states
tic
       [states_data,seed]=sim_states(X0,N,dt,type,theta,varargin{1});
disp('state simulation took:')
toc
%% Price derivatives
tic
        % Strike ranges - the range of RELATIVE strike prices to be considered at each time point
        KK = [1 1.1]'; 

        % Maturities - vector of maturities considered to create option prices at each maturity point
        TT = [0.1 0.5 1 ]'; 

        % Initialise output matrix
        deriv_data(:,1)=reshape(repmat(1:1:N,[length(TT).*length(KK),1]),[length(KK).*length(TT).*N,1]);
        aux=repmat(TT',[length(KK),1]);
        deriv_data(:,2)=repmat(aux(:),[N,1]);
        aux=repmat(KK,[1,length(TT)]); 
        deriv_data(:,3)=repmat(aux(:),[N,1]);
        deriv_data(:,6)=zeros;
        % Initialize CF coefficient matrix
%         for i=1:length(TT)
           coeff_mat=price_fft_fut(0,0,TT,0,type,theta,0,1); 
%         end

        %Pricing call options
        for i=1:N   
            %Pricing part
            for j=1:length(TT)
            aux=(price_fft_fut(states_data(i,:)', KK.*exp(states_data(i,1)), TT(j),0,type,theta,1,0,coeff_mat(:,:,j)));
            deriv_data(1+(i-1)*length(KK)*length(TT)+(j-1)*length(KK):length(KK)*length(TT)*(i-1)+j*length(KK),4:5)=[aux(:,2),exp(states_data(i,1)).*ones([length(KK),1])];
            end
        end
%         name_for_save=sprintf('deriv_data_%3.0f_%5.0f.mat',N,seed);
%         save(sprintf('Data Sim/Deriv Samples/%s',name_for_save));

% Determine and output the derivative sample structure
        M=max(deriv_data(:,1));
        index=1;
        i=1;    
        while i<=M
          TT_marks(i,1)=i;
          TT_marks(i,2)=0;
          %Pick the options in each day
          in_day=1;
          in_TT=3;
          TT_now=deriv_data(index,2);
          while index<=length(deriv_data) && deriv_data(index,1)<=i
          if(TT_now<deriv_data(index,2))
              TT_marks(i,in_TT)=in_day-1;
              in_TT=in_TT+1;
              TT_now=deriv_data(index,2);
          end
          index=index+1;
          in_day=in_day+1;
          end  
          TT_marks(i,in_TT)=in_day-1;
          i=i+1;
        end
disp('pricing options took:');
toc
end

