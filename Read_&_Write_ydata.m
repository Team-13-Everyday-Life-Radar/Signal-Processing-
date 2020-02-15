frames = 2; % chirps per frame
samples = 4; % samples per chirp
receivers = 1; % text write and read doesn't work for receivers = 2
current_line = 0;

% Made up ydata matrix
% R = rand(4,2,2);
% I = 1i*rand(4,2,2);
% ydata_write = R + I;

% % % Upload ydata matrix to text file (can only upload receiver 1 data)
% % (This algorithm doesn't work with two complex number columns)
% fileID = fopen('ydata.txt','w');
% for i = 1:frames
%     fprintf(fileID,'%f %f\n',[real(ydata_write(:,1,i)),imag(ydata_write(:,1,i))].');
% end
% fclose(fileID);

% % Read ydata matrix from text file (can only upload receiver 1)
ydata_read = zeros(samples,receivers,frames)
input = load('ydata.txt')
for h = 1:(frames)
    for g = 1:samples
        ydata_read(g,1,h) = input(g+current_line,1) + 1i*input(g+current_line,2);
    end
    current_line = current_line+g;
end
disp(ydata_read) % print ydata to command line (only contains information for receiver 1
