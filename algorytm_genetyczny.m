format compact; clc; clear; close all
%% Genetic Algorithm
global dx x_l x_u res elitist_prob initial_popu selection_prob mutation_prob
x_l = -5; % lower boundry
x_u = 15; % upper boundry
res = 10; % bit resolution of coding (2^res = length of chromosome)
dx = (x_u-x_l)/(2^res-1)

elitist_prob = 0.5
initial_popu = 20
mutation_prob = 0.5
selection_prob = 1 - elitist_prob

max_iterations = 100
min_iterations = 20

%% Goal function
gf1 = @(x) x.^2 - 10.*x + 4 - 20*sin(x);
gf1_x_min = 7.610;

figure(1)
plot (x_l:(x_u-x_l)/(2^res-1):x_u,gf1(x_l:(x_u-x_l)/(2^res-1):x_u),'k-')

%hold on

population = GenerateInitialPopulation();
i=0;
x_min = [x_u x_l];

dx_min = max(x_min)-min(x_min);
x_min_hist = [x_min];
while (i<max_iterations && dx<=dx_min)
    i=i+1;;;;;;;;;;;
     
      
    %x = DeCode(population);
    %y = gf1(x); %FittnesFunction(population,gf1) 
    %plot (x,y,'bx')
    
    plot (x_l:(x_u-x_l)/(2^res-1):x_u,gf1(x_l:(x_u-x_l)/(2^res-1):x_u),'k-','DisplayName','f(x)')
    xlabel('x')
    ylabel('f(x)')
    title(sprintf('P_{mutation}=%g, P_{elitist}=%g, P_{select}=%g, N_{population}=%g, i_{max}=%g, i_{min}=%g, \\Delta=%g',mutation_prob,elitist_prob,selection_prob,initial_popu,max_iterations,min_iterations,dx))
    hold on
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
    hold off
    if (length(x_min) < min_iterations)
        x_min = [x_min GetMinimum(population,gf1)];
    else
        x_min = [x_min(2:end) GetMinimum(population,gf1)];
    end
    dx_min = max(x_min)-min(x_min);;;;;;;;;
    x_min_hist = [x_min_hist x_min(end)];
    legend
    pause(0.01)
    input('x')
end

hold on
 x = GetMinimum(population,gf1)
 y = gf1(x); %FittnesFunction(population,gf1)
 plot (x,y,'ro','DisplayName','MINIMUM')
hold off

figure(2)
hold on
%plot(gf1(x_min_hist(3:end)),'DisplayName','MINIMUM')
plot(x_min_hist(3:end),'DisplayName','MINIMUM')
xlabel('iteration')
%ylabel('f(x_{min})')
ylabel('x_{min}')
title(sprintf('P_{mutation}=%g, P_{elitist}=%g, P_{select}=%g, N_{population}=%g, i_{max}=%g, i_{min}=%g, \\Delta=%g',mutation_prob,elitist_prob,selection_prob,initial_popu,max_iterations,min_iterations,dx))
%plot(gf1(gf1_x_min).*ones(i-1),'r--')
plot(gf1_x_min.*ones(i-1),'r--')

function PlotPopulation(population,gf,style)
    x = DeCode(population);
    y = gf(x);
    plot (x,y,style,'DisplayName',inputname(1)); 
end

function x_min = GetMinimum(p,gf)
    [~,index] = max(FittnesFunction(p,gf));
    x_min = DeCode(p(index));
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
    if(x<=x_l)
        c = 0;
    elseif (x>=x_u)
        c = 2.^res-1; 
    else
        c = floor((x-x_l)./(x_u-x_l)*(2.^res-1));
    end
end


%% Function returning decoded version of chromosome
% c - chromosome
function x = DeCode(c)
    global x_l x_u res
    x = x_l + (c./(2.^res-1)).*(x_u-x_l);
end
