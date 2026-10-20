function [hdr, adcData] = udp_receive_pulse(rx)
%UDP_RECEIVE_PULSE  Receive and reassemble one PRI sent by udp_send_pulse.m.
%   [hdr, adcData] = UDP_RECEIVE_PULSE(rx) reads datagrams from the udpport
%   object rx (MATLAB: udpport('datagram', 'LocalPort', port); Octave:
%   udpport('LocalPort', port)) until all packets of one PRI have arrived.
%   It is a reference implementation of the stream format for a receiver.
bytes = [];
received = [];
nPackets = inf;
while numel(received) < nPackets
    pkt = read_datagram(rx);
    if isempty(pkt), error('udp_receive_pulse: timeout'); end
    if typecast(pkt(1:4), 'uint32') ~= hex2dec('52414452'), continue; end   % not ours
    hdr.pulseNumber   = double(typecast(pkt(5:8),   'uint32'));
    hdr.scanNumber    = double(typecast(pkt(9:12),  'uint32'));
    hdr.antennaAz_deg = double(typecast(pkt(13:16), 'single'));
    w16 = double(typecast(pkt(17:24), 'uint16'));  % index, count, channels, stagger
    w32 = double(typecast(pkt(25:32), 'uint32'));  % samples per channel, sample offset
    hdr.time_s        = typecast(pkt(33:40), 'double');
    hdr.staggerIndex  = w16(4);
    nPackets = w16(2);
    nChannels = w16(3);
    nSamples = w32(1);
    if isempty(bytes)
        bytes = zeros(1, 2*nSamples*nChannels, 'uint8');
    end
    payload = pkt(41:end);
    bytes(2*w32(2) + (1:numel(payload))) = payload;
    received = unique([received, w16(1)]);
end
adcData = reshape(typecast(bytes, 'int16'), nSamples, nChannels);
end

function pkt = read_datagram(rx)
% One datagram as a uint8 row vector (empty on timeout).
if exist('OCTAVE_VERSION', 'builtin')
    t0 = tic;
    while rx.NumBytesAvailable == 0 && toc(t0) < rx.Timeout
        pause(0.001);
    end
    if rx.NumBytesAvailable == 0, pkt = []; return; end
    pkt = uint8(read(rx, rx.NumBytesAvailable, 'uint8'));
else
    dg = read(rx, 1, 'uint8');
    if isempty(dg), pkt = []; return; end
    pkt = uint8(dg.Data);
end
pkt = pkt(:).';
end
