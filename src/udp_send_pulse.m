function udp_send_pulse(u, hdr, ddcData)
%UDP_SEND_PULSE  Send the DDC output of one PRI as a sequence of UDP datagrams.
%   UDP_SEND_PULSE(u, hdr, ddcData) sends the complex DDC output ddcData
%   [nSamples x nChannels] of the PRI described by hdr (see pulse_header in
%   run_radar_sim.m). Every datagram = 64-byte header + up to u.payloadBytes
%   bytes of samples. Header (little-endian):
%     bytes  0-3   uint32  magic 0x52444443 ('RDDC')
%     bytes  4-7   uint32  pulse number
%     bytes  8-11  uint32  scan number
%     bytes 12-15  single  antenna azimuth [deg]
%     bytes 16-19  uint16  packet index (0-based), packet count
%     bytes 20-23  uint16  number of channels, stagger index
%     bytes 24-31  uint32  samples per channel, first sample in this packet (0-based)
%     bytes 32-39  double  PRI start time [s]
%     bytes 40-47  single  DDC sample rate [Hz], time of the first DDC sample in the PRI [s]
%     bytes 48-63  uint32  short pulse first/last DDC sample, long pulse first/last DDC sample (1-based)
%   Payload: float32 pairs (I, Q), channel after channel (column-major order of ddcData).
iq = single([real(ddcData(:)).'; imag(ddcData(:)).']);
bytes = typecast(iq(:).', 'uint8');
perPacket = 8*floor(u.payloadBytes/8);               % whole complex samples
nPackets = ceil(numel(bytes)/perPacket);
for i = 1:nPackets
    first = (i-1)*perPacket + 1;
    last  = min(i*perPacket, numel(bytes));
    head = [typecast(uint32(hex2dec('52444443')), 'uint8'), ...
            typecast(uint32([hdr.pulseNumber, hdr.scanNumber]), 'uint8'), ...
            typecast(single(hdr.antennaAz_deg), 'uint8'), ...
            typecast(uint16([i-1, nPackets, size(ddcData, 2), hdr.staggerIndex]), 'uint8'), ...
            typecast(uint32([size(ddcData, 1), (first-1)/8]), 'uint8'), ...
            typecast(double(hdr.time_s), 'uint8'), ...
            typecast(single([hdr.ddcFs_Hz, hdr.ddcT0_s]), 'uint8'), ...
            typecast(uint32([hdr.shortDdcIdx, hdr.longDdcIdx]), 'uint8')];
    packet = [head, bytes(first:last)];
    if strcmp(u.kind, 'udpport')
        write(u.obj, packet, 'uint8', u.host, u.port);
    else
        dp = javaObject('java.net.DatagramPacket', typecast(packet, 'int8'), numel(packet), u.addr, u.port);
        u.obj.send(dp);
    end
end
end
