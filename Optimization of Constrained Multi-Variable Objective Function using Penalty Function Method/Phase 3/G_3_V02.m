clc;
clear all;
%%Programming Project Phase III ME609
%
%By
%Aman Kumar
%Roll No. 214103404
%Tahir Manuel D’Mello
%Roll No. 180106055

file_input;
[extrema, function_value, results_table] = bracket_operator(10);

extrema
function_value
results_table

graphs(results_table);
                
function file_input %Strip parameters from text file 

    global problem_number;
    global variables_number;
    global epsilon;
    global x_initial;
    global epsilon_penalty;
    
    fileID = fopen('input_3.txt','r'); %Change file number here for testing convenience. File numbers correspond to question number.
    % input_1  
    % input_2  
    % input_3 
    
    ip = fscanf(fileID,'%f'); 

    problem_number = round(ip(1)); %Saves problem number as integer
    variables_number = round(ip(2)); %Saves variable number as integer
    epsilon = 10^(-ip(3)); %Saves epsilon as power of 10
    epsilon_penalty = 10^(-ip(4)); %Saves epsilon as power of 10
    
    x_initial = []; %Reinitialize set to clear previous run 
        
    for i = 1 : variables_number
        x_initial(i) = ip(4 + i); %Saves each element of initial guess as needed 
    end
    x_initial
end  

function fvalue = objfn_optimize(x_in, s_in, alpha, R)
    global feval;
    global problem_number;
    
    feval = feval + 1;
    
    x = x_in + alpha*(s_in); %Takes the input x(k), s(j) and alpha and calculates the actual x that is to be substituted into the equations
    
    if problem_number == 1
        fvalue = (x(1)-10)^3 + (x(2)-20)^3;
        fvalue = fvalue + R*constrv(x,2);
        
    elseif problem_number == 2
        fvalue = (((sin(2*pi*x(1)))^3)*(sin(2*pi*x(2))))/(x(1)^3*(x(1)+x(2)));
        fvalue = -fvalue + R*constrv(x,2); %For maximization
       
    elseif problem_number == 3
        fvalue = x(1) + x(2) + x(3);
        fvalue = fvalue + R*constrv(x,2);
        
    end
    
end %Calculate objective function using constant x_in, s_in, R and varying alpha

function out = constrv(x, type)
    global problem_number
    
    if problem_number == 1
        x1 = x(1);
        x2 = x(2);

        c(1) = min((x1-5)^2 + (x2-5)^2 - 100,0);
        c(2) = min(82.81 - (x1-6)^2 - (x2-5)^2,0);
        c(3) = min((x1-13),0);
        c(4) = min((20-x1),0);
        c(5) = min(x2,0);
        c(6) = min(4-x2,0);

    elseif problem_number == 2
        x1 = x(1);
        x2 = x(2);

        c(1) = min(-x1^2 + x2 - 1,0);
        c(2) = min(-1 + x1 - (x2-4)^2,0);
        c(3) = min(10-x1,0);
        c(4) = min(10-x2,0);
        c(5) = min(x1,0);
        c(6) = min(x2,0);

    elseif problem_number == 3
        x1 = x(1);
        x2 = x(2);
        x3 = x(3);
        x4 = x(4);
        x5 = x(5);
        x6 = x(6);
        x7 = x(7);
        x8 = x(8);

        c(1) = min(1-0.0025*(x4+x6),0);
        c(2) = min(1-0.025*(-x4+x5+x7),0);
        c(3) = min(1-0.01*(-x6+x8),0);
        c(4) = min(-(100*x1 + (-x6*x1) + (833.33252*x4) + (-83333.333)),0);
        c(5) = min(-(x2*x4 + (-x2*x7) + (-1250*x4) + (1250*x5)),0);
        c(6) = min(-(x3*x5 + (-x3*x8) + (-2500*x5) + 1250000),0);
        c(7) = min(10000-x1,0);
        c(8) = min(10000-x2,0);
        c(9) = min(10000-x3,0);
        c(10) = min(x1-100,0);
        c(11) = min(x2-1000,0);
        c(12) = min(x3-1000,0);
        c(13) = min(1000-x4,0);
        c(14) = min(1000-x5,0);
        c(15) = min(1000-x6,0);
        c(16) = min(1000-x7,0);
        c(17) = min(1000-x8,0);
        c(18) = min(x4-10,0);
        c(19) = min(x5-10,0);
        c(20) = min(x6-10,0);
        c(21) = min(x7-10,0);
        c(22) = min(x8-10,0);

    end

    if type == 1 %Check if point is feasible or infeasible
        out = logical(sum(c));
    elseif type == 2 %Returns Bracket Operator Penalty value
        out = sum(c.^2);
    elseif type == 3 %Returns constraint violation
        out = sum(c);
    end

end %Calculate penalty. Type 1 - Feasibile or not, type 2 - bracket operator value, type 3 - constraint violation

function fvalue = objfn(x, R, type)
    %global feval;
    %feval = feval + 1;
    global problem_number
    
    if problem_number == 1
        fvalue = (x(1)-10)^3 + (x(2)-20)^3;
        if type == 1
            fvalue = fvalue + R*constrv(x,2);
        end
    
    elseif problem_number == 2
        fvalue = (((sin(2*pi*x(1)))^3)*(sin(2*pi*x(2))))/(x(1)^3*(x(1)+x(2)));
        if type == 1
            fvalue = -fvalue + R*constrv(x,2);
            fvalue = -fvalue;
        end
        
    elseif problem_number == 3
        fvalue = x(1) + x(2) + x(3);
        if type == 1
            fvalue = fvalue + R*constrv(x,2);
        end
    end
    
end %Calculate objective function using x, R. Type 0 - Without penalty, Type 1 - With penalty

function out = linear_independence(s_set) 
    s_test = transpose(s_set); %Makes the row wise directions of the input into columns to calculate the rank easier
    global variables_number;
    
    if (rank(s_test) == variables_number)
        out = 1; %if rank = number of columns (no. of variables), then it is independent
    elseif(rank(s_test) < variables_number)
        out = 0; %if rank < number of columns (no. of variables), then some or all directions are dependent
    end
    
end  %Checks if any directions are linearly independent

function [x,y] = bounding_phase(x_in, s_in, R)
       
    %Delta
    delta = 0.1; %Can be changed
    k = 0;  
    
    %Initial Guess
    a0 = 1;  
    %Arbitrarily chooses a alpha intial guess for bounded phase optimization. Set to 1 for convenience.
    %Can be chosen using rand funtion too.
      
    
    %Determining the Sign of Delta
    a3 = a0 - delta;
    a4 = a0 + delta;
    
    
    %Calculates needed objfunction wrt to an alpha by sending in a fixed x_in, s_in but varying alpha. 
    %This is effectively a substitution for different alpha values
    f0 = objfn_optimize(x_in, s_in, a0, R); 
    f3 = objfn_optimize(x_in, s_in, a3, R);
    f4 = objfn_optimize(x_in, s_in, a4, R);
    
    
    if (f3 >= f0 && f0 >= f4)
        delta = delta;
    else
        delta= (-1)*delta;
    end

    %Initializing Calcuation for bounding phase method
    i = 1;
    
    a_old = a0; 
    fold = f0; 
    a_new = a_old + (2^k)*delta; %New alpha calculated
    
    fnew = objfn_optimize(x_in, s_in, a_new, R); %Objfunction for same x_in and s_in but different alpha.
    
    
    a_store = a_old; %Old alpha saved
    
    while (fold > fnew)
        
        a_store = a_old;  %Old alpha saved
        
        a_old = a_new; %New alpha made old
        fold = fnew;
               
        a_new = a_old + (2^k)*delta;  %New alpha created
        fnew = objfn_optimize(x_in, s_in, a_new, R);
              
       
        k = k + 1;
        
        i = i + 1;
                
    end
   
    if(a_store < a_new) %Output range
        x = a_store;
        y = a_new;      
    else
        x = a_new;
        y = a_store;   
    end
   
end %Gives bounding phase alpha results for a certain x_in, s_in and R

function extrema = golden_section(a, b, x_in, s_in, R)
        
    %Intitializing
    k = 1;
    anew = a;
    bnew = b;
    
    aw = (anew - a)/(b  - a); %aw = 0
    bw = (bnew - a)/(b  - a); %bw = 1
    Lw = bw - aw; %Lw = 1
    
    w1 = aw + 0.618 * Lw; %w1 and w2 calculated
    w2 = bw - 0.618 * Lw;
    
    %Normalizing and calculating objective function with respect to alpha.
    %x_in and s_in are constants here so the only variable changing is alpha.
    %Thus, the input to obj function is normalized wrt alpha
    f1 = objfn_optimize(x_in, s_in, a + w1*(b-a), R); %Normalizing equation with respect to alpha 
    f2 = objfn_optimize(x_in, s_in, a + w2*(b-a), R);%Normalizing equation with respect to alpha
    
    feval = 2;
   
    global epsilon;
    
    while Lw > epsilon
        
        %Region Elimination based on function values
        if (f1 <= f2)
           aw = w2;
           bw = bw;
           Lw = bw - aw;
           
           w1 = aw + 0.618*Lw;
           w2 = bw - 0.618*Lw;
           
           f2 = f1;
           f1 = objfn_optimize(x_in, s_in, a + w1*(b-a), R); %Normalizing equation with respect to alpha
           
        elseif(f2 < f1)
           bw = w1;
           aw = aw;
           Lw = bw - aw;
           
           w1 = aw + 0.618*Lw;
           w2 = bw - 0.618*Lw;
           
           f1 = f2;
           f2 = objfn_optimize(x_in, s_in, a + w2*(b-a), R); %Normalizing equation with respect to alpha
        end
        
        k = k + 1;
        feval = feval + 1;
    end
    
   if(f1 < f2) %Choosing needed extrema output
        extrema = a + w1*(b-a); 
    else
        extrema = a + w2*(b-a); 
   end
    
end %Gives golden section alpha results for a certain x_in, s_in and R

function extrema = powell_conjugate(x, R)

    global variables_number;
    global epsilon
    
    
    %Setting up first 'n' independent search directions using identity matrix of size n
    I = eye(variables_number);
    for i = 1 : variables_number
        s(i,:) = I(i,:);
    end
    
    k = 1;
    x(k,:) = x; %Saving initial guess to x(1,:)
    
    s_new = epsilon + 1;
    iter = 0;
    
    %Starting N^2 unidirectional search loops
    %Terminate if new direction is bigger than epsilon or the new search direction is linearly dependent
    while ( norm(s_new) > epsilon && linear_independence(s) == 1 )
                
        for j = 1 : variables_number
            [a,b] =  bounding_phase(x(k,:), s(j,:), R);
            alpha = golden_section(a, b, x(k,:), s(j,:), R);
               
            
            x(k+1,:) = x(k,:) + alpha*s(j,:); 
            %x(k+1,:) saves the next point. When j = 1, this point is taken as y(1) for conjugate direction.
            %After y(1) is saved, 'number of variables - 1' searches are done and added.
            k = k + 1;
        end
        
        [a,b] =  bounding_phase(x(k,:), s(1,:), R);
        alpha = golden_section(a, b, x(k,:), s(1,:), R); 
        %After the last s(n,:) search, another s(1,:) seach is done on the last point.
        %The result of this search gives us y(2) for conjugate direction.
        
 
        x(k+1,:) = x(k,:) + alpha*s(1,:);
        k = k + 1;
        
        %It has been 'number of variables - 1 searches' + '1 search' = 'number of variables' searches between y(2) and y(1).
        %y(2) is saved at x(k,:). 
        %Thus, y(1) is saved at x(k - variables_number,:).
        %y(2) - y(1) gives us the new conjugate direction.
        s_new = x(k,:) - x(k - variables_number,:);
        
        for j = variables_number:-1:2
            s(j,:) = s(j-1,:); %s(j,:) is made s(j-1,:)
        end
        
        s(1,:) = s_new/norm(s_new); %New direction s = y(2) - y(1) saved as s(1,:)
        s(1,:);
        
        iter = iter + 1;
    end
    
    extrema = x(k,:); %Last point identified.
    %fprintf("The minima using Powell's Conjugate method is %f. \n", round(fout,3))
    %fprintf("The minima using Powell's Conjugate method is found to be at [ ")
    %fprintf('%g  ', extrema);
    %fprintf(']\n');
    
    %fprintf("\nThe number of iterations is %d. \n", iter)
    %fprintf("The initial guess was: [ ")
    %fprintf('%g  ', x_initial);
    %fprintf(']\n');
    
end %owell Conjugate Functions for initial x and R

function [extrema, function_value, results_table] = bracket_operator(n) 
    
    global feval;
    global x_initial;
    global epsilon_penalty;
    
    c = 100; %Can be changed
    R = 0.1;
    
    %Intial Values
    feval = 0;
    results{1,1} = 0;
    results{1,2} = R; %R
    results{1,3} = x_initial; %New Point
    results{1,4} = objfn(x_initial, R, 0); %Only objective function
    results{1,5} = constrv(x_initial,2); %Penalty 
    results{1,6} = objfn(x_initial, R, 1); %Objective with penalty
    results{1,7} = constrv(x_initial,3); %Constraint violation
    results{1,8} = feval; %Number of evaluations 
    
    feval = 0;
    x_old = x_initial;
    x_next = powell_conjugate(x_old, R);
    
    results{2,1} = 1;
    results{2,2} = R; %R
    results{2,3} = x_next; %New Point
    results{2,4} = objfn(x_next, R, 0); %Only objective function
    results{2,5} = constrv(x_next,2); %Penalty 
    results{2,6} = objfn(x_next, R, 1); %Objective with penalty
    results{2,7} = constrv(x_next,3); %Constraint violation
    results{2,8} = feval; %Number of evaluations    
    
    i = 3;
    
    while i ~= -1 && i < n+1
        
        feval = 0;
        x_old = x_next;
        R = c*R;
        x_next = powell_conjugate(x_old, R);
        
        results{i,1} = i-1;
        results{i,2} = R; %R
        results{i,3} = x_next; %Point
        results{i,4} = objfn(x_next, R, 0); %Only objective function
        results{i,5} = constrv(x_next,2); %Penalty 
        results{i,6} = objfn(x_next, R, 1); %Objective with penalty
        results{i,7} = constrv(x_next,3); %Constraint violation
        results{i,8} = feval; %Number of evaluations   
        
        %fprintf("%d\n",i)
        
        %(objfn(problem_number, x_next, R, 1) - objfn(problem_number, x_old, R/c, 1)
        if abs(results{i,6} - results{i-1,6}) < epsilon_penalty
            fprintf("Penalty Function Method Terminated\n");
            break;
        end
        
        i = i + 1;
        
    end
    
    
    function_value = objfn(x_next, R, 0);
    extrema = x_next;
    results_table = cell2table(results,'VariableNames',{'Iteration','R','Extrema','Objective Function','Bracket-Operator Penalty Value','Objective Function with Penalty','Constraint Violation','Function Evaluations'});
    
end %Main bracket-operator function (n is the max amount of iterations)

function graphs(in_table)
    new_ylabel = ["Objective Function","Bracket-Operator Penalty Value", "Objective Function with Penalty", "Constraint Violation", "Function Evaluations"];

    s = stackedplot(in_table{:,1}, in_table{:,4:8});
    s.xlabel('Iterations');
    s.DisplayLabels = ["Objective Function","Bracket-Operator Penalty Value", "Objective Function with Penalty", "Constraint Violation", "Function Evaluations"];
    %s.XLimits = int[0 2];
end %Plots needed graphs