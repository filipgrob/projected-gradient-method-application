clear 
close all
clc
format long e


%Initializations
c1 = 1e-4;
rho = 0.8;
btmax = 50;
gamma = 1e-1;
tolx = 1e-12;
ftosave_fw=ones(6,3);
ftosave_grad=ones(3,1);
ftosave_c=ones(3,1);
type='0';
disp('************************************')

for d=3:5
    
    n=10^d; %dimension
    d=d-2;
    f=@(x) sum((x.^2).*(1:length(x))'); %De Jong's function
    kmax=10000; 
    tolgrad=1e-12;
    
    x0=-5*ones(n,1);
    box_mins=-5.12*ones(n,1);
    box_maxs=+5.12*ones(n,1);
    Pi_X = @(x) box_projection_sol(x, box_mins, box_maxs); %projection function
    
    
    switch type
        case "fw"
            for t=2:2:12
                gradf = @(x) findiff_grad_sol(f, x, t, type);
                %RUN THE NEWTON ON De Jong function
                disp('**** CONSTR. STEEPEST DESCENT: START De Jong function *****')
                [d,t]
                disp('                    ...')
                [xk_n, fk_n, gradfk_norm_n, deltaxk_norm_n, k_n, xseq_n, btseq_n] = ...
                    constr_steepest_desc_bcktrck_sol(x0, f, gradf, ...
                    kmax, tolgrad, c1, rho, btmax, gamma, tolx, Pi_X);
                disp('**** CONSTR. STEEPEST DESCENT: FINISHED *****')
                disp('************************************')
                ftosave_fw(t/2,d)=fk_n;
            end
        case "c"
            gradf = @(x) findiff_grad_sol(f, x, 1, type); %t=1 because it isn't used in findiff_grad
            %RUN THE NEWTON ON De Jong function
            disp('**** CONSTR. STEEPEST DESCENT: START De Jong function *****')
            d
            disp('                    ...')
            [xk_n, fk_n, gradfk_norm_n, deltaxk_norm_n, k_n, xseq_n, btseq_n] = ...
                constr_steepest_desc_bcktrck_sol(x0, f, gradf, ...
                kmax, tolgrad, c1, rho, btmax, gamma, tolx, Pi_X);
            disp('**** CONSTR. STEEPEST DESCENT: FINISHED *****')
            disp('************************************')
            ftosave_c(d)=fk_n;
        otherwise 
            gradf=@(x) 2*x.*(1:length(x))';
            %RUN THE NEWTON ON De Jong function
            disp('**** CONSTR. STEEPEST DESCENT: START De Jong function *****')
            d
            disp('                    ...')
            [xk_n, fk_n, gradfk_norm_n, deltaxk_norm_n, k_n, xseq_n, btseq_n] = ...
                constr_steepest_desc_bcktrck_sol(x0, f, gradf, ...
                kmax, tolgrad, c1, rho, btmax, gamma, tolx, Pi_X);
            disp('**** CONSTR. STEEPEST DESCENT: FINISHED *****')
            disp('************************************')
            ftosave_grad(d)=fk_n;
    end
end

%Results
switch type
    case 'fw'
        [(2:2:12)', ftosave_fw]
    case 'c'
        [(3:5)', ftosave_c] 
    otherwise
        [(3:5)', ftosave_grad]
end

%%
%Plotting results
vec1k = [4.289450618017081e-01,1.392176015056341e-17,1.392360548060494e-17,1.392305798957194e-17,1.392305253198299e-17,1.392305247740916e-17];
vec10k=[2.206331815289144e+02,7.247839916827115e-04,8.660282707176550e-02,8.659704342485960e-02,8.659699627447327e-02,8.659699580403851e-02];
vec100k=[1.258005811025823e+11,1.536864570502553e+00,3.839905308451603e+00,3.859172802575781e+00,3.859163084787475e+00,3.859162991915454e+00];
vec_c=[1.392305247685789e-17,8.659699579929085e-02,3.859162990977660e+00];
semilogy(2:2:12,vec1k,'b',2:2:12,vec10k,'r',2:2:12,vec100k,'g')
title('Plot min(f(x)) in the fw case','FontSize',12);
legend('n=1000','n=10000','n=100000');
xlabel('Value of k')
ylabel('min f(x)')

%Second plot
loglog([1000,10000,100000],vec_c,'b',[1000,10000,100000],[min(vec1k),min(vec10k),min(vec100k)],'r')
title('Plot min(f(x)) comparing fw case with c case','FontSize',12);
legend('type!=fw ', 'type=fw')
xlabel('Dimensions')
ylabel('min f(x)')