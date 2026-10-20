function ok = test_udp_stream()
%TEST_UDP_STREAM  Send one PRI over UDP (localhost) and check that it is received intact.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
port = 50123;
cfgUdp = struct('host', '127.0.0.1', 'port', port, 'payloadBytes', 8192);
if exist('OCTAVE_VERSION', 'builtin')
    pkg load instrument-control
    rx = udpport('LocalPort', port, 'Timeout', 2);
else
    rx = udpport('datagram', 'LocalPort', port, 'Timeout', 2);
end
tx = udp_sender_open(cfgUdp);

adcData = int16(randi([-3000 3000], 400, 48));            % small PRI: fits the receive socket buffer
hdr = struct('pulseNumber', 1234, 'scanNumber', 2, 'antennaAz_deg', 123.5, ...
             'staggerIndex', 3, 'time_s', 12.345);
udp_send_pulse(tx, hdr, adcData);
[hdrRx, dataRx] = udp_receive_pulse(rx);
udp_sender_close(tx);
clear rx

ok = isequal(dataRx, adcData) && hdrRx.pulseNumber == 1234 && hdrRx.scanNumber == 2 ...
     && abs(hdrRx.antennaAz_deg - 123.5) < 1e-4 && hdrRx.staggerIndex == 3 && hdrRx.time_s == 12.345;
fprintf('UDP stream test (%d bytes in %d datagrams): %s\n', 2*numel(adcData), ...
        ceil(2*numel(adcData)/cfgUdp.payloadBytes), ternary(ok, 'PASS', 'FAIL'));
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
