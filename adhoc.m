function [ diffs_map ] = adhoc( dist )

% room : 6.4, 4.5 
%  1 0 2 
%  7   3
%  6 5 4
% resolution: 0.2m 

diffs_map = zeros(45, 64);

for x = 0.1:0.1:6.4 
   for y = 0.1:0.1:4.5 
       propose_dist_lst = [
           sqrt((x-3.2)^2 + (y-0)^2);
           sqrt((x-0)^2 + (y-0)^2);
           sqrt((x-6.4)^2 + (y-0)^2);
           sqrt((x-6.4)^2 + (y-2.25)^2);
           sqrt((x-6.4)^2 + (y-4.5)^2);
           sqrt((x-3.2)^2 + (y-4.5)^2);
           sqrt((x-0)^2 + (y-4.5)^2);
           sqrt((x-0)^2 + (y-2.25)^2);
       ];
       diffs = sum(((propose_dist_lst - dist').^2)./dist'); 
       diffs_map(uint8(y*10), uint8(x*10)) = diffs;
   end
end

end

