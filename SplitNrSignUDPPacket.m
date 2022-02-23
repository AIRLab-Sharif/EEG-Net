function [Header,Data] = SplitNrSignUDPPacket(Packet,Normalized) 
    PacketSize = size(Packet,1);
    if (PacketSize < 3)
        Header = {};
        Data = [];        
    else
        SamplesCount = Packet(1) + 256 * Packet(2);
        Header = cell(1,65); 
        Data = zeros(SamplesCount,65);
        Index = 3;
        i = 1;
        while (Index<=PacketSize)
            HeaderLength = Packet(Index);
            Header{i} = char(Packet(Index+1 : Index+HeaderLength))'; 
            if Normalized % in uV
                Index = Index + HeaderLength + 1;
                for j = 1:SamplesCount
                    Data(j,i) = typecast(uint8(Packet(Index:Index+3)), 'single');
                    Index = Index + 4;
                end               
            else % Raw Data
                Data(:,i) = Packet(Index+HeaderLength+1 : 2 : ...
                     Index+HeaderLength+1+2*(SamplesCount-1)) ...
                     + 256 * Packet(Index+HeaderLength+2 : 2 : ... 
                     Index+HeaderLength+2+2*(SamplesCount-1));
                Index = Index+HeaderLength+2+2*(SamplesCount-1)+1;
            end            
            i=i+1;
        end
        if ~Normalized
            Data = Data - 32768;
        end
    end
end