Normalized = true; % false/true, depending on software online-setting
HeaderIndex = [1, 2]; % Indices of the headers/channels to be plotted
SampleRate = 500; % Hz/Per second (depending on the device)
PlotLength = 5; % In seconds; Minimum: Buffer_Size/Sample_Rate.
PlotRefreshRate = 15; % Plot refreshing rate in Hz (<= 15 Hz).

HeaderNumber = length(HeaderIndex); % Number of channels to be plotted
BufferLength = PlotLength*SampleRate; % data length to be plotted 
BufferData = zeros(BufferLength,HeaderNumber); 

t = [0:1/SampleRate:(BufferLength-1)/SampleRate]'; % Time (s)
y = zeros(BufferLength,1); % Output signal

ServerPort = 12220;
ClientPort = 12221;

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
tic; toc_old=toc;
while(ByteCount)
    % Reading new data
    [Packet,ByteCount] = fread(u,1);
    [Header,Data] = SplitNrSignUDPPacket(Packet,Normalized);
    size(Data)

    if isempty(Data)
        disp('No data is available.');
    else
        if JustStarted 
            %%% Plotting initialization (2) %%%               
            for hd=1:HeaderNumber
                temp = min(8,ceil(HeaderNumber/2));
                subplot(temp,ceil(HeaderNumber/temp),hd);
                plt(hd) = plot(t,y);
                %ylim([-1000 1000]);
                title(Header{HeaderIndex(hd)});
                plt(hd).XDataSource = 't';
                plt(hd).YDataSource = 'y';
            end
            JustStarted = false;
        end
        for hd = 1 : HeaderNumber
            % Storing new data in buffer
            BufferData(:,hd)=[BufferData(size(Data,1)+1:end,hd); ...
                              Data(:,HeaderIndex(hd))];
            % Plotting new data   
            y=BufferData(:,hd);
            refreshdata(plt(hd),'caller');
        end

        if (toc-toc_old > 1/PlotRefreshRate)% Drawing the updated data
            drawnow;
            toc_old=toc;
        end            
    end
end
fclose(u);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
