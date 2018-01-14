% Distributed representation of fuzzy rules and its application to pattern
% classification

function Main


clc, clear

h = 0;
rate_per_each_iterations=[];

   for t=1 : 10
     
%----------------generate train data------------
    
    generate_random_value_data_train = rand(100,2);
    x1 = generate_random_value_data_train(:,1);
    x2 = generate_random_value_data_train(:,2);
    problem_value = (-1/4) * sin(2* pi * x1) + x2 - (0.5);
%   problem_value = (-1/3) * sin(2* pi * x1) + x2 - (0.5);
%   problem_value = (-1/3) * sin(2* pi * x1 - (1/2) * pi ) + x2 - (0.5);
%   problem_value = (-1) * abs((-2) * x1 + 1) + x2;  
%   problem_value = (x1 + x2 - 1).*(x2 - x1);
%   problem_value = ((-(x1 - 0.5).^2)/(0.4).^2) - (((x2 - (0.5)).^2)/(0.3).^2) + 1;
%   problem_value = ((-(x1 - 0.5).^2)/(0.15).^2) + (((x2 - (0.5)).^2)/(0.2).^2) + 1;
    
    
     for i= 1 : size (problem_value)
       if   problem_value(i) >= 0
               problem_value(i) = 1;
       else    problem_value(i) = 2;
           
       end
     end
     
% ----------------plot train data----------------     
%     figure, hold
%     for i=1:size(problem_value)
%         z = 'r.';
%         if( problem_value(i) == 2) 
%             z = 'b.';
%         end
%         plot(x1(i), x2(i), z);
%     end



%------------------variables definition--------------
   data_train_matric = [x1,x2,problem_value];
   rules_region_g1 = [];
   rules_class_g1  = 0;
   rules_cf_g1     = [];
   rules_k_g1      = [];
   rules_region_g2 = [];
   rules_class_g2  = 0;
   rules_cf_g2     = [];
   rules_k_g2      = [];
   rules_region_unclass = [];
   rules_unclass   = 0;
   rules_k_unclass = [];
   q = 0;
   g1 = data_train_matric(data_train_matric(:,3) == 1, :);
   g2 = data_train_matric(data_train_matric(:,3) == 2, :);
   L_Maximum = 20;
   
   
   
%---------------create fuzzy rules---------------------- 
    for k = 2 : L_Maximum
%        fprintf('\nfor k = %i rules are :\n', k);  
      for i = 1 : k
         for j = 1 : k
            q = q + 1;
            beta_class_g1 = sum(triangular(g1(:,1), i, k) .* triangular(g1(:,2), j, k));
            beta_class_g2 = sum(triangular(g2(:,1), i, k) .* triangular(g2(:,2), j, k));
%           beta_class_g1 = sum(min(triangular(g1(:,1), i, k), triangular(g1(:,2), j, k)));
%           beta_class_g2 = sum(min(triangular(g2(:,1), i, k), triangular(g2(:,2), j, k)));
%           beta_class_g1 = sum(trapezoid(g1(:,1), i, k) .* trapezoid(g1(:,2), j, k));
%           beta_class_g2 = sum(trapezoid(g2(:,1), i, k) .* trapezoid(g2(:,2), j, k));
            cf = abs(beta_class_g1 - beta_class_g2) / (beta_class_g1 + beta_class_g2);
            if beta_class_g1 > beta_class_g2
%                fprintf('if "i" is %i and "j" is %i then class is G%i and CF = %f\n', i, j, 1, cf);
                rules_region_g1(end+1,:) = [i j];
                rules_class_g1       = rules_class_g1 + 1;
                rules_cf_g1(end + 1) = cf;
                rules_k_g1(end+1)    = k;
            elseif beta_class_g2 > beta_class_g1   
%                fprintf('if "i" is %i and "j" is %i then class is G%i and CF = %f\n', i, j, 2, cf);
                rules_region_g2(end+1,:) = [i j];
                rules_class_g2       = rules_class_g2 + 1;
                rules_cf_g2(end + 1) = cf;
                rules_k_g2(end+1)    = k;
            else
%                fprintf('if "i" is %i and "j" is %i then unclass\n', i, j);
                rules_region_unclass(end+1,:)  = [i j];
                rules_unclass     = rules_unclass + 1;
%               rules_cf(end + 1) = cf;
                rules_k_unclass(end+1)=k;
            end
        end
      end
    end
%    fprintf('number of generated rules are : %i \n',q);
   
   
%-----------------generate test data--------------------   
    
    generate_random_value_for_test_set = rand(100,2);
    y1 = generate_random_value_for_test_set(:,1);
    y2 = generate_random_value_for_test_set(:,2);
    test_data_problem = (-1/4) * sin(2* pi * y1) + y2 - (0.5);
%   test_data_problem = (-1/3) * sin(2* pi * y1) + y2 - (0.5);
%   test_data_problem = (-1/3) * sin(2* pi * y1 - (1/2) * pi ) + y2 - (0.5);
%   test_data_problem = (-1) * abs((-2) * y1 + 1) + y2;  
%   test_data_problem = (y1 + y2 - 1).*(y2 - y1);
%   test_data_problem = ((-(y1 - 0.5).^2)/(0.4).^2) - (((y2 - (0.5)).^2)/(0.3).^2) + 1;
%   test_data_problem = ((-(y1 - 0.5).^2)/(0.15).^2) + (((y2 - (0.5)).^2)/(0.2).^2) + 1;
    
     for i= 1 : size (test_data_problem)
       if test_data_problem(i) >= 0;
              test_data_problem(i) = 1;
       else   test_data_problem(i) = 2;
       end
     end
     
     t= test_data_problem;
    
     
%--------------plot test data---------------------------     
%     figure, hold
%     for i=1:size(tast_date_problem)
%         z = 'r.';
%         if( tast_date_problem(i) == 2) 
%             z = 'b.';
%         end
%         plot(y1(i), y2(i), z);
%     end
%     

%-------------fuzzy inference for test data pattern classification-------
   alpha_g1 = 0; 
   for i=1 : rules_class_g1
         z = triangular(generate_random_value_for_test_set(:, 1), rules_region_g1(i, 1), rules_k_g1(i)) .* triangular(generate_random_value_for_test_set(:, 2), rules_region_g1(i, 2), rules_k_g1(i)) .* rules_cf_g1(i);
%        z = min(triangular(generate_random_value_for_test_set(:, 1), rules_region_g1(i, 1), rules_k_g1(i)) , triangular(generate_random_value_for_test_set(:, 2), rules_region_g1(i, 2), rules_k_g1(i))) .* rules_cf_g1(i);
%        z = triangular(generate_random_value_for_test_set(:, 1), rules_region_g1(i, 1), rules_k_g1(i)) .* triangular(generate_random_value_for_test_set(:, 2), rules_region_g1(i, 2), rules_k_g1(i));
%        z = trapezoid (generate_random_value_for_test_set(:, 1), rules_region_g1(i, 1), rules_k_g1(i)) .* trapezoid (generate_random_value_for_test_set(:, 2), rules_region_g1(i, 2), rules_k_g1(i)) .* rules_cf_g1(i);
%        z = trapezoid (generate_random_value_for_test_set(:, 1), rules_region_g1(i, 1), rules_k_g1(i)) .* trapezoid (generate_random_value_for_test_set(:, 2), rules_region_g1(i, 2), rules_k_g1(i));
         alpha_g1 = max(alpha_g1, z);
   end
   
   alpha_g2 = 0;
   for i=1 : rules_class_g2
         z = triangular(generate_random_value_for_test_set(:, 1), rules_region_g2(i, 1), rules_k_g2(i)) .*  triangular(generate_random_value_for_test_set(:, 2), rules_region_g2(i, 2), rules_k_g2(i)) .* rules_cf_g2(i);
%        z = min(triangular(generate_random_value_for_test_set(:, 1), rules_region_g2(i, 1), rules_k_g2(i)) ,  triangular(generate_random_value_for_test_set(:, 2), rules_region_g2(i, 2), rules_k_g2(i))) .* rules_cf_g2(i);
%        z = triangular(generate_random_value_for_test_set(:, 1), rules_region_g2(i, 1), rules_k_g2(i)) .*  triangular(generate_random_value_for_test_set(:, 2), rules_region_g2(i, 2), rules_k_g2(i));
%        z = trapezoid (generate_random_value_for_test_set(:, 1), rules_region_g2(i, 1), rules_k_g2(i)) .*  trapezoid (generate_random_value_for_test_set(:, 2), rules_region_g2(i, 2), rules_k_g2(i)) .* rules_cf_g2(i);
%        z = trapezoid (generate_random_value_for_test_set(:, 1), rules_region_g2(i, 1), rules_k_g2(i)) .*  trapezoid (generate_random_value_for_test_set(:, 2), rules_region_g2(i, 2), rules_k_g2(i));
         alpha_g2 = max(alpha_g2, z);
   end
   
   c = 0;
   results = zeros(size(generate_random_value_for_test_set,1), 2);
   for i=1 : 100
       

       if    (alpha_g1(i) > alpha_g2(i)), c = 1;
       elseif(alpha_g2(i) > alpha_g1(i)), c = 2;
       end
       results(i,:) = [c , abs(alpha_g2(i) - alpha_g1(i))];
   
   end
   
%   results
    
    
%------------------rate computation------------------
    u = results(:,1);
    r = zeros();
    r = t - u ;
    v = 0;
     for i = 1:size(generate_random_value_for_test_set)
         if  r(i) == 0;
             v = v + 1;
         end
     end

rate_per_each_iterations(end+1)= v;
h = h + v;




%---------------- plot result---------------------
%     figure, hold
%     for i=1:size(results,1)
%         z = 'r.';
%         if( results(i, 1) == 2), z = 'b.';
%         elseif(results(i, 2) == 0), z = 'k.';
%         end
%         plot(generate_random_value_for_test_set(i, 1), generate_random_value_for_test_set(i, 2), z);
%     end



   end
  
  rate_per_each_iterations
  total_rate = h/10
 
 end




%-------------------membership functions-------------------------
function membership_value = triangular(x,i,k)
a = (i-1)/(k-1);
b = 1/(k-1);
membership_value = max(1- abs(x-a)/b,0);
end

function membership_value = trapezoid(x,i,k)
a = (i-1)/(k-1);
b = 1/(k-1);
membership_value = max(min(2 - 2 * abs(x-a) / b, 1));
end


