%% New Algorithm IRM with increased sampling in each radius
%load('D:\Research\Sem 9\Model Assessment\Examples\A linear Biochemical pathway\xd.mat')
%% Sampling parameters

mu=[1 10]; % Optimal parameter set

delta = 0:0.01:0.3; % Radius of n-ball around the optimal parameters

%% Generating parameter samples around the paramter of interest

N0=100000;% initial sample size
alpha=50000;
m1=[];
m1(2)=N0;

for i=2:length(delta)
    
    z=1;
    gamma_temp=[];
 
    temp_draw=[];
    
    Para_sample=[];
    
    n=2; % number of dimentions
    
    r = delta(i); % radius of the n-sphere
    
    if i ==2
        
        m=N0;
        
    else
        m1(i)= m1(i-1) + alpha*(delta(i)/delta(i-1))^n;% number of samples
        
        m=round(m1(i));
        
    end
    
    
    X = randn(m,n);
    
    s2 = sum(X.^2,2);
    
    X = X.*repmat(r*(rand(m,1).^(1/n))./sqrt(s2),1,n); % Projecting the points on to the sphere
    
    Para_sample = X+mu;
    % toc
    %% Solving ODE
    
    
    for l=1:length(Para_sample)
        
        temp_para=Para_sample(l,:);
        
        del_val(i)=delta(i-1);
        
        n_temp(l)=norm(temp_para-mu);
        
        if  n_temp(l) > delta(i-1)
            
            %% Simulating ODE
            temp_draw(z,:)=temp_para;
            
            tspan=0:0.25:30;% time span
            
            p=temp_para;
            
             
              %x0=[0.3617 0.9137 1.3934];
            
            %opts = odeset('AbsTol',1e-3);
            
            %[t,x]=ode23s(@(t,x)pk1(t,x,p),tspan,x0);% ode solver
            x=exp(-p(1)*tspan)+ exp(-p(2)*tspan);
            
              %gamma_temp(z)= norm(data(:,1)-x(:,1));
            
            gamma_temp(z)= norm(x-data); % Measure of closeness
            
            z=z+1;
            
        end
    end
    
    
%    [gamma_min(i), I(i)]=min(gamma_temp(i,:));% Computing gamma-min
%    [gamma_max(i), I1(i)]=max(gamma_temp(i,:));% Computing gamma-max
    
    [gamma_min(i),I(i)] =min(gamma_temp);% Computing gamma-min
    [gamma_max(i),I1(i)]=max(gamma_temp);% Computing gamma-max
    
    %
    Min_vecs(i,:)=temp_draw(I(i),:);% Extracting Sloppy Vector
    Min_norms(i)=norm(mu-Min_vecs(i,:));
    Max_vecs(i,:)=temp_draw(I1(i),:);% Extracting Stiff Vector
    Max_norms(i)=norm(mu-Max_vecs(i,:));
end
%save(Algo_inc_sample)
