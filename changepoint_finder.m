function [cp_index,slopes] = changepoint_finder(radii,distance)

n = length(radii);

if n < 10
    cp_index = [];
    slopes = [];
else
    
%One changepoint
best1_error1 = Inf;
best1_error2 = Inf;

for point = 3:n-2

data1 = 1:point;
data2 = point:length(distance);

oldWarningState = warning;
warning('off','all');
p1 = polyfit(distance(data1),radii(data1),1);

p2 = polyfit(distance(data2),radii(data2),1);

seg1 = polyval(p1,distance(data1));
seg2 = polyval(p2,distance(data2));

cost_line = sum((radii(data1)-seg1').^2)+...
                  sum((radii(data2)-seg2').^2);
cost_mean = sum((radii(data1) -mean(radii(data1))).^2) +...
            sum((radii(data2) -mean(radii(data2))).^2);

    if (cost_line < best1_error1) && (cost_mean < best1_error2)
        cp1 = point;
        best1_error1 = cost_line;
        best1_error2 = cost_mean;
        seg1_1 = seg1;
        seg1_2 = seg2;

        slopes1 = [p1(1) p2(1)];
                    
         
        params1_seg1 = [mean(radii(data1)) var(radii(data1))];
        params1_seg2 = [mean(radii(data2)) var(radii(data2))];

        loglik1_seg1 = -0.5 * log(2 * pi * params1_seg1(2)) - 0.5 *...
                               ((radii(data1) - params1_seg1(1)).^2) / params1_seg1(2);
        loglik1_seg2 = -0.5 * log(2 * pi * params1_seg2(2)) - 0.5 *...
                               ((radii(data2) - params1_seg2(1)).^2) / params1_seg2(2);
        
        log_likes1 = [loglik1_seg1; loglik1_seg2];                  
        loglik1 = sum(log_likes1); 
        
    end
            
end

% 2  changepoints
best2_error1 = Inf;
best2_error2 = Inf;

for j=3:n-2
    for k=j+3:n-2
        data1 = 1:j;
        data2 = j:k;
        data3 = k:length(radii);
                
        p1 = polyfit(distance(data1),radii(data1),1);
        p2 = polyfit(distance(data2),radii(data2),1);
        p3 = polyfit(distance(data3),radii(data3),1);

        seg1 = polyval(p1,distance(data1));
        seg2 = polyval(p2,distance(data2));
        seg3 = polyval(p3,distance(data3));
        
      
       cost2_line = sum((radii(data1)-seg1').^2)+...
                    sum((radii(data2)-seg2').^2)+...
                    sum((radii(data3)-seg3').^2);
       cost2_mean = sum((radii(data1) -mean(radii(data1))).^2) +...
                    sum((radii(data2) -mean(radii(data2))).^2) +...
                    sum((radii(data3) -mean(radii(data3))).^2);


              

          if (cost2_line < best2_error1) && (cost2_mean < best2_error2)
            cp2 = [j,k];
            best2_error1 = cost2_line;
            best2_error2 = cost2_mean;
            
            seg2_1 = seg1;    
            seg2_2 = seg2;
            seg2_3 = seg3;
            
            slopes2 = [p1(1) p2(1) p3(1)];
                        
            
            params2_seg1 = [mean(radii(data1)) var(radii(data1))];
            params2_seg2 = [mean(radii(data2)) var(radii(data2))];
            params2_seg3 = [mean(radii(data3)) var(radii(data3))];
            
            loglik2_seg1 = -0.5 * log(2 * pi * params2_seg1(2)) - 0.5 *...
                           ((radii(data1) - params2_seg1(1)).^2) / params2_seg1(2);
            loglik2_seg2 = -0.5 * log(2 * pi * params2_seg2(2)) - 0.5 *...
                           ((radii(data2) - params2_seg2(1)).^2) / params2_seg2(2);
            loglik2_seg3 = -0.5 * log(2 * pi * params2_seg3(2)) - 0.5 *...
                           ((radii(data3) - params2_seg3(1)).^2) / params2_seg3(2); 
                       
            log_likes2 = [loglik2_seg1; loglik2_seg2; loglik2_seg3];                  
            loglik2 = sum(log_likes2);     
            
            
            
         end
    end                   
end

%3 changepoints
best3_error1 =Inf;
best3_error2 =Inf;

for j=3:n-3
    for k=j+3:n-2
        for l=k+3:n-1
        data1 = 1:j;
        data2 = j:k;
        data3 = k:l;
        data4 = l:length(radii);
        
        p1 = polyfit(distance(data1),radii(data1),1);
        p2 = polyfit(distance(data2),radii(data2),1);
        p3 = polyfit(distance(data3),radii(data3),1);
        p4 = polyfit(distance(data4),radii(data4),1);

        seg1 = polyval(p1,distance(data1));
        seg2 = polyval(p2,distance(data2));
        seg3 = polyval(p3,distance(data3));
        seg4 = polyval(p4,distance(data4));              
        
        cost3_line = sum((radii(data1)-seg1').^2)+...
                     sum((radii(data2)-seg2').^2)+...
                     sum((radii(data3)-seg3').^2)+...
                     sum((radii(data4)-seg4').^2);
                 
        cost3_mean = sum((radii(data1) -mean(radii(data1))).^2) +...
                     sum((radii(data2) -mean(radii(data2))).^2) +...
                     sum((radii(data3) -mean(radii(data3))).^2) +...
                     sum((radii(data4) -mean(radii(data4))).^2);

        if  (cost3_line < best3_error1) && (cost3_mean < best3_error2)
             
            cp3 = [j,k,l];
            best3_error1 = cost3_line;
            best3_error2 = cost3_mean;
            
            seg3_1 = seg1;
            seg3_2 = seg2;
            seg3_3 = seg3;
            seg3_4 = seg4;
            
            
            slopes3 = [p1(1) p2(1) p3(1) p4(1)];

            params3_seg1 = [mean(radii(data1)) var(radii(data1))];
            params3_seg2 = [mean(radii(data2)) var(radii(data2))];
            params3_seg3 = [mean(radii(data3)) var(radii(data3))];
            params3_seg4 = [mean(radii(data4)) var(radii(data4))];
            
            loglik3_seg1 = -0.5 * log(2 * pi * params3_seg1(2)) - 0.5 *...
                           ((radii(data1) - params3_seg1(1)).^2) / params3_seg1(2);
            loglik3_seg2 = -0.5 * log(2 * pi * params3_seg2(2)) - 0.5 *...
                           ((radii(data2) - params3_seg2(1)).^2) / params3_seg2(2);
            loglik3_seg3 = -0.5 * log(2 * pi * params3_seg3(2)) - 0.5 *...
                           ((radii(data3) - params3_seg3(1)).^2) / params3_seg3(2);
            loglik3_seg4 = -0.5 * log(2 * pi * params3_seg4(2)) - 0.5 *...
                           ((radii(data4) - params3_seg4(1)).^2) / params3_seg4(2);           
            
            log_likes3 = [loglik3_seg1; loglik3_seg2; loglik3_seg3; loglik3_seg4];                  
            loglik3 = sum(log_likes3);  

        end
 
       end
    end
end


loglikelihood = [loglik1, loglik2, loglik3];
penalty = log(length(radii));

for k = 1:3 %up to 3 changepoints
    bic(k) = k * penalty - 2 * loglikelihood(k);
    aic(k) = 2*k - 2*loglikelihood(k);
    bic_fake(k) =  2*(k-1)*penalty- (loglikelihood(k));
end

num_cp_bic = find(bic==min(bic));
num_cp_aic = find(aic==min(aic));
num_cp_fake = find(bic_fake==min(bic_fake));

cp_num = num_cp_fake;


if cp_num == 3
    cp_index = cp3;
    slopes = slopes3;
elseif cp_num == 2
    cp_index = cp2;
    slopes = slopes2;
else 
    cp_index = cp1;
    slopes = slopes1;
end
    
end