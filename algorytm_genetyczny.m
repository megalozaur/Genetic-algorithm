format compact; clc; clear; close all;
%% Genetic Algorithm
global dx x_l x_u res elitist_prob initial_popu selection_prob mutation_prob debug
x_l = [-5,0,0]; % lower boundry
x_u = [15,3,2]; % x_u - upper boundry
res = length(x_l)*3; % bit resolution of coding (2^res = length of chromosome)
dx = min(x_u-x_l)/(2^(res/length(x_l))-1)
elitist_prob = 0.2
initial_popu = 40
mutation_prob = 0.2
selection_prob = 1 - elitist_prob
 debug=10;

max_iterations = 1000
min_iterations = 100

%% Goal function
gf1 = @(x) x(:,1).^2 - 10.*x(:,1) +  20.*x(:,1).*sin(x(:,1)) + 4 - 20*sin(x(:,2)) + 10*x(:,3);
gf1_min = -224.9291
PlotFunction(gf1)
population = GenerateInitialPopulation();
i=0;
y_min = [gf1(x_u) gf1(x_l)];

dy_min = max(y_min)-min(y_min);
y_min_hist = [y_min];
while (i<max_iterations && dx<=dy_min)
    i=i+1;;;;;;;;;;;
     
    PlotFunction(gf1)
    
    PlotPopulation(population,gf1,'bx')
    
    survivors = ElitistSelection(population,gf1);;;;;;;;;;
    PlotPopulation(survivors,gf1,'bo')
    
    parents = ParentsSelection(population,gf1);;;;;;;;;;;;;
    allparents = [parents(1,:) parents(2,:)];
    PlotPopulation(allparents,gf1,'r^')
  
    newborns = OnePointCrossover(parents);;;;;;;;;;;;;
    PlotPopulation(newborns,gf1,'gs')
    mutants = Mutate(newborns);;;;;;;;;;;;;;;
    
    %mutants = new_borns
    PlotPopulation(mutants,gf1,'mp')
    population = [survivors mutants];
    %population = mutants;
    
    [~,~,y_min_i] = GetMinimum(population,gf1);
    if (length(y_min) < min_iterations)
        y_min = [y_min y_min_i];
    else
        y_min = [y_min(2:end) y_min_i];
    end
    dy_min = max(y_min)-min(y_min);;;;;;;;;
    y_min_hist = [y_min_hist y_min(end)];
    
end
i
[x,c,y] = GetMinimum(population,gf1)
if (debug>0)
    hold on
    plot (c,y,'ro','DisplayName','MINIMUM')
    hold off
end

figure(2)
hold on
%plot(gf1(x_min_hist(3:end)),'DisplayName','MINIMUM')
plot(y_min_hist(3:end),'DisplayName','MINIMUM')
xlabel('iteration')
ylabel('f(x_{min})')

title(sprintf('P_{mutation}=%g, P_{elitist}=%g, P_{select}=%g, N_{population}=%g, i_{max}=%g, i_{min}=%g, \\Delta=%g',mutation_prob,elitist_prob,selection_prob,initial_popu,max_iterations,min_iterations,dx))
plot(gf1_min.*ones(i-1),'r--')

function PlotFunction(gf)
    global dx res x_l debug elitist_prob initial_popu mutation_prob selection_prob max_iterations min_iterations
    if (debug>0)
        figure(1)
        title(sprintf('P_{mutation}=%g, P_{elitist}=%g, P_{select}=%g, N_{population}=%g, i_{max}=%g, i_{min}=%g, \\Delta=%g',mutation_prob,elitist_prob,selection_prob,initial_popu,max_iterations,min_iterations,dx))
        hold on
        c = 0:1:(2^res-1);
        y = gf(DeCode(c));
        plot (c,y,'k-','DisplayName','f(x)')
        xlabel('coded x')
        ylabel('f(x)')
        legend
    end
end

function PlotPopulation(population,gf,style)
    global debug
    if (debug>0)
        hold on
        x = population;
        y = gf(DeCode(x));
        plot (x,y,style,'DisplayName',inputname(1)); 
        hold off
    end
end

function [x_min,c,y] = GetMinimum(p,gf)
    [y,index] = min(-FittnesFunction(p,gf));
    c = p(index);
    x_min = DeCode(c);
end

function mutants = Mutate(p)
    global res mutation_prob
    % pattern = randi(2.^res,1,length(p))-1; % mutation_prob = 0.5
    pattern = sum((rand(1,res) < mutation_prob) .* 2.^(res-1:-1:0));
    mutants = bitxor(p,pattern);
end

function new_borns = OnePointCrossover(parents)
    global res
    % parents(1,:) - ojciec
    % pparents(2,:) - matka
    cut_points = randi(res-1,1,size(parents,2));
    pattern_1 = 2.^cut_points - 1;
    pattern_2 = 2^res - pattern_1 - 1;
    new_borns = bitand(parents(1,:),pattern_1) + bitand(parents(2,:),pattern_2);
    
end

function parents = ParentsSelection(p,gf)
    global selection_prob
    n = floor(selection_prob*length(p)); % number of pairs of chrom. to select
    if (n<=1)
        n=1;
    end
    
    [FF,indexes] = sort(FittnesFunction(p,gf),'descend');
    if(max(FF)==min(FF)) 
        parents=[p(1:n);p(1:n)];
    else
        series_of_ranges = (FF - min(FF)) ./ (max(FF)-min(FF));
        % https://la.mathworks.com/matlabcentral/answers/51897-generate-random-numbers-given-distribution-histogram
        D = cumsum(series_of_ranges)/sum(series_of_ranges);
        f = @(r) find(r<=D,1,'first');
        parents = p( arrayfun(f,rand(n,2)) )';
    end
end

function p = ElitistSelection(p,gf)
    global elitist_prob
    m = elitist_prob*length(p); % number of chrom. to select
    if (m<1)
        m=1;
    end
    [~,indexes] = sort(FittnesFunction(p,gf),'descend');
    p = p(indexes);
    p = p(1:m);
end

function p = GenerateInitialPopulation()
    global res initial_popu
    p = randi(2.^res-1,1,initial_popu);
end

%% Fittness Function
% gf - goal function of x
% c - chromosome
function r = FittnesFunction(c,gf)
    r = - gf(DeCode(c));
end

%% Function returning coded version of value
% x - value from range x_l : x_u
function c = Code(x)
    global x_l x_u res
    x_norm = (x-x_l)./(x_u-x_l);
    x_int = floor(x_norm*(2.^(res/length(x_l))-1));
    weights = 2.^(0:length(x)*(res/length(x_l))-1);
    c = sum(x_int(:)*weights)
end

%% Function returning decoded version of chromosome
% c - chromosome
function x = DeCode(c)
    global x_l x_u res
    xyz_bin = de2bi(c,res);
    x_bin = reshape(xyz_bin',res/length(x_l),length(x_l)*length(c));
    x_int=bi2de(x_bin');
    xyz_int = reshape(x_int,length(x_l),length(c))';
    x_norm = xyz_int/(2.^(res/length(x_l))-1);
    x = (x_l + x_norm.*(x_u-x_l));
end
