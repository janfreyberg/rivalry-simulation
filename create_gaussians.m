%% Random Variables (location, slope)
stimsize = 500;
transition_frames = 90;
n = 8;
[x, y] = meshgrid(linspace(-stimsize, stimsize, stimsize));
while 1
    mu_x = randi([0.1*stimsize, 0.9*stimsize], [1, 10]);
    mu_y = randi([0.1*stimsize, 0.9*stimsize], [1, 10]);
    small_distance = stimsize;
    for i = 1:n
        
        other_ns = 1:10;
        other_ns(other_ns==i) = [];
        
        for ii = other_ns
            if sqrt( (mu_x(i) - mu_x(ii))^2 + (mu_y(i) - mu_y(ii))^2 ) < small_distance
                small_distance = sqrt( (mu_x(i) - mu_x(ii))^2 + (mu_y(i) - mu_y(ii))^2 );
                sminx = mu_x([i, ii]);
                sminy = mu_y([i, ii]);
            end 
        end
        
    end
    if small_distance > 40
        break
    end
end

max_stdev = stimsize/3; % this is how large the gaussians are
max_amplitude = 100;

gauss_onset = rand(1, n)*stimsize/15; % the blobs start randomly in the first fifth of the transition period

gauss_slope = 0.5 + rand(1, n);

stdev_x = linspace(0, 100, transition_frames);

amp = zeros(transition_frames,n);
for i=1:n
    amp(:, i) = gauss_slope(1,i).*(stdev_x-gauss_onset(1,i));
end

amp(amp>max_amplitude) = max_amplitude;
amp(amp<0)=0;
amp = amp/100;
figure;
subplot(1, 3, 1);
plot(amp);
subplot(1, 3, 2);
scatter(mu_x, mu_y);



for ii = 1:transition_frames
disp(ii)
gaussians(ii).red = zeros(500,500);
for i = 1:n
% tmp = (3*amp(ii,i)/stimsize).*exp(-(((x-mu_x(i))/(stimsize/3)).^2)-(((y-mu_y(i))/(stimsize/3)).^2));
% gaussians(ii).red(tmp>gaussians(ii).red) = tmp(tmp>gaussians(ii).red);
gaussians(ii).red(:,:) = gaussians(ii).red(:,:) + amp(ii,i).*exp(-(((x-mu_x(i))/(stimsize/3)).^2)-(((y-mu_y(i))/(stimsize/3)).^2));
end
gaussians(ii).red(:,:) = gaussians(ii).red(:,:)*255;
gaussians(ii).red(gaussians(ii).red>255) = 255;
end


gaussians(transition_frames).red = 255*ones(500,500);


for ii = 1:transition_frames
disp(ii)
vlm = trapz(y(:,1),trapz(x(1,:),gaussians(ii).red(:,:),2),1);
desired_vlm = (2*stimsize)^2*(ii*255/transition_frames);
graph(ii, 1) = vlm;
graph(ii, 2) = desired_vlm;
while abs(vlm-desired_vlm) > desired_vlm*0.000001
gaussians(ii).red(:,:) = ((gaussians(ii).red(:,:))/vlm) * desired_vlm;
gaussians(ii).red(gaussians(ii).red>255) = 255;
vlm = trapz(y(:,1),trapz(x(1,:),gaussians(ii).red(:,:),2),1);
end
end
subplot(1, 3, 3);
plot(graph);
