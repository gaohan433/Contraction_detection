function [] = show_stat(metric_array)
    x1 = metric_array(1,:);
    fprintf('median x1 = %.4f\n', median(x1));
    fprintf(' IQR x1 = %.4f\n', iqr(x1))
    for i = 1:size(metric_array,1)-1
        x2 = metric_array(i+1,:);
        disp(i)
        fprintf('median x2 = %.4f\n', median(x2));
        fprintf(' IQR x2 = %.4f\n', iqr(x2))
        [p,h,stats] = signrank(x1,x2);
        
        diffs = x1 - x2;
  
        diffs(diffs == 0) = [];
        N = numel(diffs);
        
        % Compute Z statistic
        W = stats.signedrank;
        Z = (W - (N*(N+1))/4) / sqrt(N*(N+1)*(2*N+1)/24);
        
        % Compute effect size
        r = Z / sqrt(N);
        
        fprintf('p = %.4f, Z = %.4f, r = %.4f\n', p, Z, abs(r));
        disp('------')
    end
end


