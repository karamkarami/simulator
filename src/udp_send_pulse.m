function udp_send_pulse(u, hdr, adcData)
%UDP_SEND_PULSE  Send the raw ADC data of one PRI as a sequence of UDP datagrams.
%   Every datagram = 40-byte header + up to u.payloadBytes bytes of int16 samples.
%   Header (little-endian):
%     uint32 magic 0x52414452 ('RADR') | uint32 pulse number | uint32 scan number
%     single antenna azimuth [deg]     | uint16 packet index (0-based) | uint16 packet count
%     uint16 number of channels        | uint16 stagger index
%     uint32 samples per channel       | uint32 offset of the first sample in this packet
%     double PRI start time [s]
%   Samples are sent channel after channel (column-major order of adcData).
bytes = typecast(adcData(:), 'uint8');
nPerPacket = 2*floor(u.payloadBytes/2);
nPackets = ceil(numel(bytes)/nPerPacket);
for i = 1:nPackets
    first = (i-1)*nPerPacket + 1;
    last  = min(i*nPerPacket, numel(bytes));
    head = [typecast(uint32(hex2dec('52414452')), 'uint8'), ...
            typecast(uint32(hdr.pulseNumber), 'uint8'), ...
            typecast(uint32(hdr.scanNumber), 'uint8'), ...
            typecast(single(hdr.antennaAz_deg), 'uint8'), ...
            typecast(uint16([i-1, nPackets, size(adcData, 2), hdr.staggerIndex]), 'uint8'), ...
            typecast(uint32([size(adcData, 1), (first-1)/2]), 'uint8'), ...
            typecast(double(hdr.time_s), 'uint8')];
    packet = [head, bytes(first:last).'];
    if strcmp(u.kind, 'udpport')
        write(u.obj, packet, 'uint8', u.host, u.port);
    else
        dp = javaObject('java.net.DatagramPacket', typecast(packet, 'int8'), numel(packet), u.addr, u.port);
        u.obj.send(dp);
    end
end
end
