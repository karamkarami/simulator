function u = udp_sender_open(cfgUdp)
%UDP_SENDER_OPEN  Open a UDP sender for streaming the raw ADC data.
%   Uses udpport (MATLAB R2020b+, or the Octave instrument-control package)
%   and falls back to a Java DatagramSocket.
u = cfgUdp;
if exist('OCTAVE_VERSION', 'builtin')
    pkg load instrument-control
end
try
    u.obj = udpport('datagram');            % MATLAB
    u.kind = 'udpport';
catch
    try
        u.obj = udpport();                  % Octave instrument-control
        u.kind = 'udpport';
    catch
        u.obj  = javaObject('java.net.DatagramSocket');
        u.addr = javaMethod('getByName', 'java.net.InetAddress', u.host);
        u.kind = 'java';
    end
end
end
