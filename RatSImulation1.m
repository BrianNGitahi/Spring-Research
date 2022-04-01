clear all;

% Rat foraging in a 1D box simulatio

%1. Define the circular arena with vectors that are <= norm 1

m = 1000;           % number of attempted points
z = 1000;           % max size of coords array
coords = randn(z,2);% array of x & y coords


% 2.1 Generate random vectors in this 2d space and put in the array if
% they fall within the circle of rad 1 centered around the origin

for i = 1:m
    
      xi = randn(1,1); 
      yi = randn(1,1); 
      
      % condition to make sure that the points are within the bounds
      while sqrt(xi*xi + yi*yi) > 1
          xi = randn(1,1); 
          yi = randn(1,1); 
      end
      
      pt = [xi, yi];  
      coords(i,:) = pt;
   
end


% 2.2 Create a sequence of points to follow while respecting the conditions set

n = 500;
sequence = randn(n,2);

% choose a random vector to start from 

rindex = randi(m);      % random index
ab = coords(rindex,:);  % random pt 
sequence(1,:) = ab;
used_pts = zeros(m, 1); % keeping track of used points
used_pts(rindex,:) = 1; % update tracker to reflect that 1st pt has been used
u = 1;                  % keeps track of index that's been used in the coords array : usedindex
s = 0.1;                % Initialize s parameter


% get subsequent points
for i = 1:n

    st = rand*s;
    currpt = sequence(i,:);

    % updated version of coords without the used indices
    coords_unused = coords(used_pts ~= 1,:);

    % finding dist between current point and all others that haven't been used in the array of points
    dist = sqrt(sum((coords_unused - currpt).^2,2));

    % Sorted version of the value array (the diff between st and dist)
    % with the indices also rearranged
    value = abs(st - dist);
    [val2, ind] = sort(value);
   
   
    l = length(ind); %% added by Manuel
   
    % iterating thru sorted version of value
   for j = 1:l 
     
       %condition that makes sure we move forward  while accounting for first point case
       if i>1
            cond1 = sum((coords_unused(ind(j),:) - currpt).^2) < sum((coords_unused(ind(j),:) - sequence(i-1,:)).^2);
       else
            cond1 = true;
       end

       flag = cond1;
         
           
  % once flag == true, then we?ve found our desired next point
       if flag
           dpt = coords_unused(ind(j),:);  % desired point
           break;                   % then break from the for loop here, right? i.e once we've found
                                    % our desired pt
       end
   end

if j == l
    dpt = coords_unused(ind(1),:);
end
   
sequence(i+1,:) = dpt;
% actual_idx = find(dpt==coords);
actual_idx = find ( sum( ( coords - dpt).^2 ) == 0);
used_pts(actual_idx,:) = 1;

disp(i)

end


figure(1)
plot(coords(:,1),coords(:,2),'.')
hold on;
plot(sequence(:,1), sequence(:,2), '-')
color = linspace(1,10,length(sequence));
scatter(sequence(:,1), sequence(:,2), [], color, 'o')


