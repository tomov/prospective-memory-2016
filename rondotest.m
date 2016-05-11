poolobj = gcp;
poolobj = gcp;

poolobj.NumWorkers
fprintf('num of workers = %d\n', poolobj.NumWorkers);

a = 1:100;

parfor i=a
    i
end


init_par = [0.35    0.3, ...     % focal, low emph 
            0.6     0.4, ...     % focal, high emph 
            0.8     0.75, ...    % nonfocal, low emph 
            0.9     0.83, ...    % nonfocal, high emph
            0.4]; % gamma * 10^-3

error = fitme(init_par);

error

fprintf('error = %f\n', error);
