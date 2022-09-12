clear all;
clc;
%%Programming Project Phase II ME609
%EXTREME POINT CALCULATION OF A FUNCTION
%By
%Tahir Manuel D’Mello
%Roll No. 180106055
%Aman Kumar
%Roll No. 214103404

%Change file number in file_input for testing convenience @Aman

global problem_number;
global variables_number;
global epsilon;
global x_initial;
global lower_bound;
global upper_bound;
global accuracy_measure;

file_input;

answer = powell_conjugate;
restart_function(answer);
                
function file_input %Strip parameters from text file 

    global problem_number;
    global variables_number;
    global epsilon;
    global x_initial;
    global lower_bound;
    global upper_bound;
    
    fileID = fopen('input_5.txt','r'); %Change file number here for testing convenience
    ip = fscanf(fileID,'%f'); 

    problem_number = round(ip(1)) %Saves problem number as integer
    variables_number = round(ip(2)); %Saves variable number as integer
    epsilon = 10^(-ip(3)); %Saves epsilon as power of 10
    
    x_initial = []; %Reinitialize set to clear previous run 
    lower_bound = [];
    upper_bound = [];
    
    for i = 1 : variables_number
        lower_bound(i) = ip(4); %Saves each element of bound as needed
        upper_bound(i) = ip(5); %Saves each element of bound as needed
        x_initial(i) = ip(5 + i); %Saves each element of initial guess as needed
    end
    
end  

function restart_function(answer)
    global x_initial;
    
    decision = input("\nDo you wish to restart the method with the current output minima as new initial guess? \nInput 1 for 'Yes' and 0 for 'No' \n");
    answer_2 = answer;
    
    while (decision == 1) %
        x_initial = answer_2;
        fprintf('\nNew initial guess is: ['); 
        fprintf('%g  ', x_initial);
        fprintf(']\n');
        
        answer_2 = powell_conjugate;
        
        decision = input("\nDo you wish to restart the method again with the current output minima as new initial guess? \nInput 1 for 'Yes' and 0 for 'No' \n");
    end
        
end %Restarts function with inital guess as last output. Iterates till '0' input.

function fvalue = objfn(problem_number, variables_number, x_in, s_in, alpha)
    
    x_in; 
    s_in;
    x = x_in + alpha*(s_in); %Takes the input x(k), s(j) and alpha and calculates the actual x that is substituted into the equations
    fvalue = 0;
    
    if (problem_number == 0)
        %fprintf("Himmenblau Function - Test. /n")
        
        fvalue = (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2; %Use only the first two elements of x
        fvalue;
      
        
    elseif (problem_number == 1)
        %fprintf("Sum of square functions. /n")
              
        for i = 1 : variables_number
            fvalue = fvalue + i*x(i)^2;
        end
        
        fvalue;
       
        
    elseif (problem_number == 2)
         %fprintf("Rosnebrock Function. /n")
                
        for i = 1 : variables_number - 1
            fvalue = fvalue + 100*(x(i+1) - x(i)^2)^2 + (x(i) - 1)^2;
        end
        fvalue;

        
    elseif (problem_number == 3)
        %fprintf("Dixon Price Function. /n")
            
        fvalue_1 = (x(1) - 1)^2;
        
        for i = 2 : variables_number
            fvalue = fvalue + i*(2*x(i)^2 - x(i-1))^2;
        end
        
        fvalue = fvalue_1 + fvalue;
        fvalue;

        
    elseif (problem_number==4)
        %fprintf("Trid Function. /n")
        
        fvalue_1 = 0;
        fvalue_2 = 0;
        
        for i = 1 : variables_number
            fvalue_1 = fvalue_1 + (x(i)-1)^2;
        end
        
        for i = 2 : variables_number
            fvalue_2 = fvalue_2 + x(i)*x(i-1);
        end
        
        fvalue = fvalue_1 - fvalue_2;
        fvalue;
      
        
    elseif (problem_number==5)
        %fprintf("Zakharov Function. /n")
        
        fvalue_1 = 0;
        fvalue_2 = 0;
        fvalue_3 = 0;
        
        for i = 1 : variables_number 
            fvalue_1 = fvalue_1 + x(i)^2;
            fvalue_2 = fvalue_2 + 0.5*i*x(i);
            fvalue_3 = fvalue_3 + 0.5*i*x(i);
        end
        
        fvalue_2 = (fvalue_2)^2;
        fvalue_3 = (fvalue_3)^4;
        fvalue = fvalue_1 + fvalue_2 + fvalue_3;
        
        fvalue;
        
    end
    
end %Calculate objective function using constant x_in, s_in and varying alpha

function out = linear_independence(s_set) 
    s_test = transpose(s_set); %Makes the row wise directions of the input into columns to calculate the rank easier
    global variables_number;
    
    if (rank(s_test) == variables_number)
        out = 1; %if rank = number of columns (no. of variables), then it is independent
    elseif(rank(s_test) < variables_number)
        out = 0; %if rank < number of columns (no. of variables), then some or all directions are dependent
    end
    
end  %Checks if any directions are linearly independent

function [x,y] = bounding_phase(x_in, s_in)
    
    global problem_number;
    global variables_number;
    global lower_bound;
    global upper_bound;
    
    %Delta
    n = 100; %Number of steps, can be changed
    delta = (upper_bound(1)-lower_bound(1))/n; %Delta calculated based on range, arbitrary but gives structure
    k = 0;  
    
    %Initial Guess
    a0 = 1; %rand(1)*(upper_bound(1) - lower_bound(1)) + lower_bound(1); 
    %Arbitrarily chooses a alpha intial guess for bounded phase optimization. Set to 1 for convenience.
    %Can be chosen using rand funtion too.
      
    
    %Determining the Sign of Delta
    a3 = a0 - delta;
    a4 = a0 + delta;
    feval = 0;
    
    %Calculates needed objfunction wrt to an alpha by sending in a fixed x_in, s_in but varying alpha. 
    %This is effectively a substitution for different alpha values
    f0 = objfn(problem_number, variables_number, x_in, s_in, a0); 
    f3 = objfn(problem_number, variables_number, x_in, s_in, a3);
    f4 = objfn(problem_number, variables_number, x_in, s_in, a4);
    feval = feval + 3;
    
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
    
    fnew = objfn(problem_number, variables_number, x_in, s_in, a_new); %Objfunction for same x_in and s_in but different alpha.
    
    feval = feval + 1;
    a_store = a_old; %Old alpha saved
    
    while (fold > fnew)
        
        a_store = a_old;  %Old alpha saved
        
        a_old = a_new; %New alpha made old
        fold = fnew;
               
        a_new = a_old + (2^k)*delta;  %New alpha created
        fnew = objfn(problem_number, variables_number, x_in, s_in, a_new);
              
        feval = feval + 1;
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
   
end %Gives bounding phase alpha results for a certain x_in and s_in

function extrema = golden_section(a, b, x_in, s_in)
    global problem_number;
    global variables_number;
    
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
    f1 = objfn(problem_number, variables_number, x_in, s_in, a + w1*(b-a)); %Normalizing equation with respect to alpha
    f2 = objfn(problem_number, variables_number, x_in, s_in, a + w2*(b-a));%Normalizing equation with respect to alpha
    
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
           f1 = objfn(problem_number, variables_number, x_in, s_in, a + w1*(b-a)); %Normalizing equation with respect to alpha
           
        elseif(f2 < f1)
           bw = w1;
           aw = aw;
           Lw = bw - aw;
           
           w1 = aw + 0.618*Lw;
           w2 = bw - 0.618*Lw;
           
           f1 = f2;
           f2 = objfn(problem_number, variables_number, x_in, s_in, a + w2*(b-a)); %Normalizing equation with respect to alpha
        end
        
        k = k + 1;
        feval = feval + 1;
    end
    
   if(f1 < f2) %Choosing needed extrema output
        extrema = a + w1*(b-a); 
    else
        extrema = a + w2*(b-a); 
   end
    
end %Gives golden section alpha results for a certain x_in and s_in

function extrema = powell_conjugate

    global variables_number;
    global x_initial;
    global epsilon
    global problem_number;
    
    %Setting up first 'n' independent search directions using identity matrix of size n
    I = eye(variables_number);
    for i = 1 : variables_number
        s(i,:) = I(i,:);
    end
    
    k = 1;
    x(k,:) = x_initial; %Saving initial guess to x(1,:)
    
    s_new = epsilon + 1;
    iter = 0;
    
    %Starting N^2 unidirectional search loops
    %Terminate is new direction is bigger than epsilon or the new search direction is linearly dependent
    while ( norm(s_new) > epsilon && linear_independence(s) == 1 )
                
        for j = 1 : variables_number
            [a,b] =  bounding_phase(x(k,:), s(j,:));
            alpha = golden_section(a, b, x(k,:), s(j,:));
               
            
            x(k+1,:) = x(k,:) + alpha*s(j,:); 
            %x(k+1,:) saves the next point. When j = 1, this point is taken as y(1).
            %After y(1) is saved, 'number of variables - 1' searches are done and added.
            k = k + 1;
        end
        
        [a,b] =  bounding_phase(x(k,:), s(1,:));
        alpha = golden_section(a, b, x(k,:), s(1,:)); 
        %After the last s(n,:) search, another s(1,:) seach is done on the last point.
        %The result of this search gives us y(2)
        
 
        x(k+1,:) = x(k,:) + alpha*s(1,:);
        k = k + 1;
        
        %It has been 'number of variables - 1' + '1 search' = 'number of variables' searches between y(2) and y(1).
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
    
    s
    x
     
    fout = objfn(problem_number, variables_number, x(k-1,:), s(1,:), alpha); %Minima identified with the last point found.
    round(fout,3);
    
    extrema = round(x(k,:),3); %Last point identified.
    fprintf("The minima using Powell's Conjugate method is %f. \n", round(fout,3))
    fprintf("The minima using Powell's Conjugate method is found to be at [ ")
    fprintf('%g  ', extrema);
    fprintf(']\n');
    
    fprintf("\nThe number of iterations is %d. \n", iter)
    fprintf("The initial guess was: [ ")
    fprintf('%g  ', x_initial);
    fprintf(']\n');
    
    fprintf("\nFor Aman: Norm of the output minima point is %f. \nI would have to change too much code to get accuracy. So just minus this with theoretical minima norm to get accuracy.\n", norm(extrema));
    
end %Main Powell Conjugate Functions.
