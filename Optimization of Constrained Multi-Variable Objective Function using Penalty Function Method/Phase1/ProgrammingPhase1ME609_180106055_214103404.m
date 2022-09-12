%%Programming Project Phase I ME609
%EXTREME POINT CALCULATION OF A FUNCTION
%By
%Tahir Manuel D’Mello
%Roll No. 180106055
%Aman Kumar
%Roll No. 214103404

%Uncomment needed function in objfn and proceed.

a=input('Enter a = lower limit of x: ');
b=input('Enter b = upper limit of x: ');
if (b<=a)
    disp('Please enter reasonable values (b>a) and reinitialize.');
    return;
end

global epsilon;
epsilon = input('Enter epsilon = small error value: ');

global t;
t = input("Please specify the type of function: 'max' or 'min'? (Give the input in lowercase only): ",'s');
if ( (strcmpi(t,'max') == 0) && (strcmpi(t,'min') == 0))
    fprintf("Please specify the type of function: 'max' or 'min' and reinitialize.\n");    
    return;
end

[x,y] = bounding_phase(a,b);
final_extrema = golden_section(x,y);



function fvalue = objfn(x)
    
    %Uncomment needed function and proceed.
    %fvalue = (2*x-5)^4-(x^2-1)^3;     %Max (-10,0) 
    %fvalue = 8+x^3-2*x-2*exp(x);        %Max (-2,1)
    fvalue = 4*x*sin(x);                %Max (0.5,pi)
    %fvalue = 2*(x-3)^2+exp(0.5*x^2);    %Min (-2,3)
    %fvalue = x^2-10*exp(0.1*x);        %Min (-6,6)
    %fvalue = 20*sin(x)-15*x^2;          %Max (-4,4)
    
    global t;
    if (t == 'max')
        fvalue = (-1)*fvalue;
    end
end

function [x,y] = bounding_phase(a,b)
    
    fprintf('\nBounding Phase Method\n');
    fprintf('--------------------------\n');
    
    %Initial Guess
    x0 = input('\nEnter an initial guess: ');
    if (x0<=a || x0>=b)
        disp('Please enter reasonable guess between upper and lower limit and reinitialize.');
        return;
    end
    
    %Delta
    n = input("Enter the number of steps: ");
    if (n<0)
        disp('Please enter n>0 and reinitialize.');
        return;
    end
    delta = (b-a)/n;
    k = 0;  
    
    %Determining the Sign of Delta
    x3 = x0 - delta;
    x4 = x0 + delta;
    feval = 0;

    f0 = objfn(x0);
    f3 = objfn(x3);
    f4 = objfn(x4);
    feval = feval + 3;
    
    if (f3 >= f0 && f0 >= f4)
        delta = delta;
    elseif (f3 <= f0 && f0 <= f4)
        delta= (-1)*delta;
    else
        fprintf('\nReinitialize with different initial guess.\n');
        return;
    end
    fprintf("\nDelta = %f \n",delta);
    
    %Initializing Calcuation
    i = 1;
    
    xold = x0; 
    fold = f0; 
    xnew = xold + (2^k)*delta;
    fnew = objfn(xnew);
    
    A(1) = fnew; %Stores values for the graph
    bounding_out = fopen('bounding_phase_iterations.out', 'w'); % Output file
    fprintf(bounding_out, '#It\tk\tx(k)\tx(k+1)\tf(x(k))\tf(x(k+1))\n');
    
    feval = feval + 1;
    
    while (fold > fnew && xnew <= b && xnew >= a)
        
        fprintf(bounding_out, '%d\t%d\t%4.4f\t%4.4f\t%4.4f\t%4.4f\n',i,k,xold,xnew,fold,fnew); %Stores iteration wise values
        
        xstore = xold;  
        
        xold = xnew;
        fold = fnew;
               
        xnew = xold + (2^k)*delta; 
        fnew = objfn(xnew);
              
        feval = feval + 1;
        k = k + 1;
        
        i = i + 1;
        A(i) = fnew;
                
    end

    fclose(bounding_out);
    
    if(xstore > b || xstore < a || xnew > b || xnew < a) %Edge case
        fprintf('\nReinitialize with different initial guess if needed, values crossing bounds. Extrema could be at boundaries.\n');
    end
    
    if(xstore < xnew)
        fprintf('\nThe extrema lies between %2.3f & %2.3f\n',xstore,xnew);
        x = xstore;
        y = xnew;      
    else
        fprintf('\nThe extrema lies between %2.3f & %2.3f\n',xnew,xstore);
        x = xnew;
        y = xstore;   
    end

    fprintf('\nTotal number of function evaluations are %d\n',feval);
    fprintf('Total iterations are %d',i-1);
    fprintf('\nThe extrema interval after applying the Bounding Phase bracketing method is (%2.3f, %2.3f)\n',x,y);
    
    %Plotting graph between Iterations (i) & function values corresponding to ith iteration.
    figure(1);
    s=[1:i];
    global t;
    if (t == 'max')
        A = (-1)*A;
    end
    plot(s,A);
    xlabel('Iterations');
    ylabel('Function Value');
    title('Bounding Phase - Function Value vs Iterations');
    
end

function extrema = golden_section(a,b)
    fprintf('\nGolden Section Method\n');
    fprintf('--------------------------\n');
    
    %Intitializing
    k = 1;
    anew = a;
    bnew = b;
    
    aw = (anew - a)/(b  - a); %aw = 0
    bw = (bnew - a)/(b  - a); %bw = 1
    Lw = bw - aw; %Lw = 1
    
    w1 = aw + 0.618 * Lw;
    w2 = bw - 0.618 * Lw;
    
    %Normalizing and calculating objective function
    f1 = objfn(a + w1*(b-a)); %Since w = (x-a)/(b-a) => x = a + w*(b-a)
    f2 = objfn(a + w2*(b-a));
    
    B(1) = f2; %Stores values for the graph
    feval = 2;
    
    golden_out = fopen('golden_section_iterations.out', 'w'); % Output file
    fprintf(golden_out, '#It\tw(1)\tw(2)\tf(w(1))\tf(w(2))\n');
    
    global epsilon;
    while Lw > epsilon
        
        fprintf(golden_out, '%d\t%4.4f\t%4.4f\t%4.4f\t%4.4f\n',k,w1,w2,f1,f2); %Stores iteration wise values
        
        %Region Elimination
        if (f1 < f2)
           aw = w2;
           bw = bw;
           Lw = bw - aw;
           
           w1 = aw + 0.618*Lw;
           w2 = bw - 0.618*Lw;
           
           f2 = f1;
           f1 = objfn(a + w1*(b-a));
           
        elseif(f2 < f1)
           bw = w1;
           aw = aw;
           Lw = bw - aw;
           
           w1 = aw + 0.618*Lw;
           w2 = bw - 0.618*Lw;
           
           f1 = f2;
           f2 = objfn(a + w2*(b-a));
        end
        
        k = k + 1;
        B(k) = f2;
        feval = feval + 1;
    end
    
    fclose(golden_out);
    
    fprintf('\nTotal number of function evaluations are %d\n',feval);
    fprintf('Total iterations are %d\n',k-1);
    
    if(f1 < f2)
        extrema = a + w1*(b-a);
        fprintf('\nThe extrema from the Golden Section Method is at %2.3f.\n',extrema);     
    else
        extrema = a + w2*(b-a);
        fprintf('\nThe extrema from the Golden Section Method is at %2.3f.\n',extrema);   
    end
      
    figure(2);
    r=[1:k];
    global t;
    if (t == 'max')
        B = (-1)*B;
    end
    plot(r,B);
    xlabel('Iterations');
    ylabel('Function Value');
    title('Golden Section - Function Value vs Iterations');
    
end
