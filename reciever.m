clc
clear
close all

load Electrodes.mat;

% configs

plt_chs1 = [2, 7, 14, 23, 32, 41, 48, 53, 56, 55, 51, 45, 37, 28, 19, 11, 5];
plt_chs2 = [1, 3, 8, 15, 24, 33, 42, 49, 54, 55, 51, 45, 37, 28, 19, 11, 5];

ch = 55;
fc = 40;
w_filter = 1.5;
order = 2;

Normalized = 1; 
HeaderIndex = 1:56;
Fs = 500;
PlotLength = 2.5; 
PlotRefreshRate = 15;

ServerPort = 12220;
ClientPort = 12221;

HeaderNumber = length(HeaderIndex); 
BufferLength = PlotLength*Fs;
BufferData = zeros(BufferLength,HeaderNumber); 

elec_names = {'FP1';'FP2';'AF7';'AF3';...
    'AFZ';'AF4';'AF8';'F7';'F5';'F1';'FZ';...
    'F2';'F6';'F8';'FT7';'FC5';'FC3';'FC1';...
    'FCZ';'FC2';'FC4';'FC6';'FT8';'T7';'C5';...
    'C3';'C1';'CZ';'C2';'C4';'C6';'T8';'TP7';...
    'CP5';'CP3';'CP1';'CPZ';'CP2';'CP4';'CP6';...
    'TP8';'P7';'P5';'P1';'PZ';'P2';'P6';'P8';...
    'PO7';'PO3';'POZ';'PO4';'PO8';'O1';'OZ';'O2'};

indexes = [];
elocsX = [];
elocsY = [];
elabels = [];

for i=1:numel(elec_names)
    fun = @(x) strcmp(Electrodes(x).labels, elec_names{i});
    tf2 = cell2mat(arrayfun(fun, 1:56, 'UniformOutput',false));
    [~,col] = find(tf2);
    indexes = [indexes; col];
    elocsX = [elocsX; Electrodes(col).X];
    elocsY = [elocsY; Electrodes(col).Y];
end

elocsY = -1*elocsY;
elabels = [elabels, ''];

[b,a] = butter(order, [fc-w_filter, fc+w_filter]/Fs*2);

t = [0:1/Fs:(BufferLength-1)/Fs]';
y = zeros(BufferLength,length(HeaderIndex));


% setting up UDP connection
if(~isempty(instrfindall))
    fclose(instrfindall);
end

u = udp('192.168.1.102','RemotePort',ServerPort, ...
    'Localport',ClientPort, 'ByteOrder','bigEndian');
set(u,'InputBufferSize',50*65535);
set(u,'InputDatagramPacketSize',65535);
set(u,'Timeout',60);
fopen(u);
ByteCount = 1;
JustStarted = true;

%%%% Main %%%%
time = 0;
while(ByteCount)
    
    % Reading new data
    [Packet,ByteCount] = fread(u,1);
    try
        [Header,Data] = SplitNrSignUDPPacket(Packet,Normalized);
        time = time+size(Data, 1)/Fs;
    catch
        ByteCount = 1;
    end
    if numel(Header)<length(HeaderIndex)
        continue
    end
    Header(9) = [];
    Data(:, 9) = [];
    Header(18) = [];
    Data(:, 18) = [];
    Header(25) = [];
    Data(:, 25) = [];
    Header(27) = [];
    Data(:, 27) = [];
    Header(end) = [];
    Data(:, end) = [];
    Header(end) = [];
    Data(:, end) = [];
    Header(end) = [];
    Data(:, end) = [];
    Header(end) = [];
    Data(:, end) = [];
    Header(end) = [];
    Data(:, end) = [];

    if isempty(Data)
        disp('No data is available.');
    else
        if JustStarted 
            %%% Plotting initialization (2) %%%               

            topo_mat = zeros(9);
            topo_mat(1, [4,6]) = y(end, 1:2);
            topo_mat(2, 3:7) = y(end, 3:7);
            topo_mat(3, [1:2, 4:6, 8:9]) = y(end, 8:14);
            topo_mat(4, 1:9) = y(end, 15:23);
            topo_mat(5, 1:9) = y(end, 24:32);
            topo_mat(6, 1:9) = y(end, 33:41);
            topo_mat(7, [1:2, 4:6, 8:9]) = y(end, 42:48);
            topo_mat(8, 3:7) = y(end, 49:53);
            topo_mat(9, 4:6) = y(end, 54:56);
            topo_mat = flip(topo_mat, 1);
            
            figure
            plt_power = imagesc(topo_mat);
            caxis([-40 0])
            colormap jet
            
            figure
            plt_wave = imagesc(topo_mat);
            caxis([-1 1])
            colormap hot

            figure
            data_ch = y(:,ch); 
            t = (1:1:length(data_ch))/Fs+time;
            plt = plot(data_ch);
            plt.XDataSource = 't';
            plt.YDataSource = 'data_ch';
            ylim([-1000, 1000])
            
            figure
            data_fft = fft(y);
            L = length(data_fft);
            f = Fs*(0:(L/2))/L;
            P2 = abs(data_fft(:, ch)/L);
            P1 = P2(1:L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);
            plt_fft = plot(f, P1);
            plt_fft.XDataSource = 'f';
            plt_fft.YDataSource = 'P1';
 
            figure
            hold on
            plts_ch1_wave = [];
            for i = 1:length(plt_chs1)
                ch_no = plt_chs1(i);
                plts_ch1_wave = [plts_ch1_wave, plot(t, (i-1)*100+y(:, ch_no))];
            end
            hold off

            figure
            hold on
            plts_ch2_wave = [];
            for i = 1:length(plt_chs2)
                ch_no = plt_chs2(i);
                plts_ch2_wave = [plts_ch2_wave, plot(t, (i-1)*100+y(:, ch_no))];
            end
            hold off

            JustStarted = false;
        end
        
        BufferData =[BufferData(size(Data,1)+1:end,:); ...
                          Data(:,HeaderIndex(:))];
        
        % Plotting in time domain
        y = BufferData;
        t = (1:1:length(data_ch))/Fs+time;
        data_ch = y(:,ch);
        refreshdata(plt,'caller');

        % Plotting in frequency domain
        data_fft = fft(y);
        f = Fs*(0:(L/2))/L;
        P2 = abs(data_fft(:, ch)/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        refreshdata(plt_fft,'caller');

        % Plotting wave
        y_wave = filter(b, a, y);
        y_wave = cos(angle(hilbert(y_wave)));
        for i = size(Data, 1):-5:1
            sample_no = BufferLength-i; 
            topo_mat(1, [4,6]) = y_wave(sample_no, 1:2);
            topo_mat(2, 3:7) = y_wave(sample_no, 3:7);
            topo_mat(3, [1:2, 4:6, 8:9]) = y_wave(sample_no, 8:14);
            topo_mat(4, 1:9) = y_wave(sample_no, 15:23);
            topo_mat(5, 1:9) = y_wave(sample_no, 24:32);
            topo_mat(6, 1:9) = y_wave(sample_no, 33:41);
            topo_mat(7, [1:2, 4:6, 8:9]) = y_wave(sample_no, 42:48);
            topo_mat(8, 3:7) = y_wave(sample_no, 49:53);
            topo_mat(9, 4:6) = y_wave(sample_no, 54:56);
            set(plt_wave, 'CData', topo_mat);
        end
        
        for i = 1:length(plts_ch1_wave)
            set(plts_ch1_wave(i), 'XData', t);
            set(plts_ch1_wave(i), 'YData', (i-1)*100+y(:,plt_chs1(i)));
        end
        
        for i = 1:length(plts_ch2_wave)
            set(plts_ch2_wave(i), 'XData', t);
            set(plts_ch2_wave(i), 'YData', (i-1)*100+y(:,plt_chs2(i)));
        end
        
        % Plotting power
        PSD = abs(data_fft).^2;
        PSD = PSD(1:L/2+1, :)/(Fs*L);
        PSD(2:end-1, :) = 2*PSD(2:end-1, :);
        PSD = 10*log10(PSD);
        [~, col] = find(f==fc);

        y = PSD(col, :);
        topo_mat(1, [4,6]) = y(:, 1:2);
        topo_mat(2, 3:7) = y(:, 3:7);
        topo_mat(3, [1:2, 4:6, 8:9]) = y(:, 8:14);
        topo_mat(4, 1:9) = y(:, 15:23);
        topo_mat(5, 1:9) = y(:, 24:32);
        topo_mat(6, 1:9) = y(:, 33:41);
        topo_mat(7, [1:2, 4:6, 8:9]) = y(:, 42:48);
        topo_mat(8, 3:7) = y(:, 49:53);
        topo_mat(9, 4:6) = y(:, 54:56);
        set(plt_power, 'CData', topo_mat);
        
        if (toc-toc_old > 1/PlotRefreshRate)% Drawing the updated data
            drawnow;
            toc_old=toc;
        end
        
    end
end
fclose(u);