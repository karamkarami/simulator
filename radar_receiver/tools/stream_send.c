/*
 * stream_send.c - sends a recorded DDC stream (radar_receiver --record FILE)
 * over UDP at the real PRI rate, to test the receiver in real time without
 * the MATLAB simulator.
 *
 *   radar_stream_send FILE [--host 127.0.0.1] [--port 50000] [--speed 1.0] [--loops 1]
 *
 * The datagrams of one PRI are sent together; the PRIs are paced by the PRI
 * start times in their headers divided by --speed (0 = as fast as possible).
 * With --loops > 1 the pulse numbers of later loops are continued, so the
 * receiver sees one long stream.
 */
#if defined(_WIN32)
#  ifndef _CRT_SECURE_NO_WARNINGS
#    define _CRT_SECURE_NO_WARNINGS
#  endif
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <winsock2.h>
#  include <ws2tcpip.h>
typedef SOCKET sock_t;
#else
#  ifndef _POSIX_C_SOURCE
#    define _POSIX_C_SOURCE 200809L
#  endif
#  include <arpa/inet.h>
#  include <netinet/in.h>
#  include <sys/socket.h>
#  include <time.h>
#  include <unistd.h>
typedef int sock_t;
#endif

#include "radar_rx.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void sleep_until(double t)
{
    const double dt = t - now_seconds();
    if (dt <= 0) return;
#ifdef _WIN32
    Sleep((DWORD)(dt * 1e3));
#else
    struct timespec ts = { (time_t)dt, (long)((dt - (double)(time_t)dt) * 1e9) };
    nanosleep(&ts, NULL);
#endif
}

static uint32_t get_u32(const unsigned char *p) { return p[0] | p[1] << 8 | p[2] << 16 | (uint32_t)p[3] << 24; }
static void put_u32(unsigned char *p, uint32_t v) { p[0] = v & 255; p[1] = v >> 8 & 255; p[2] = v >> 16 & 255; p[3] = v >> 24; }

int main(int argc, char **argv)
{
    const char *file = NULL, *host = "127.0.0.1";
    int port = 50000, loops = 1;
    double speed = 1.0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--host") && i + 1 < argc) host = argv[++i];
        else if (!strcmp(argv[i], "--port") && i + 1 < argc) port = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--speed") && i + 1 < argc) speed = atof(argv[++i]);
        else if (!strcmp(argv[i], "--loops") && i + 1 < argc) loops = atoi(argv[++i]);
        else if (argv[i][0] != '-' && !file) file = argv[i];
        else { fprintf(stderr, "usage: radar_stream_send FILE [--host H] [--port N] [--speed X] [--loops N]\n"); return 2; }
    }
    if (!file) { fprintf(stderr, "usage: radar_stream_send FILE [--host H] [--port N] [--speed X] [--loops N]\n"); return 2; }
#ifdef _WIN32
    WSADATA wsa;
    WSAStartup(MAKEWORD(2, 2), &wsa);
#endif
    sock_t s = socket(AF_INET, SOCK_DGRAM, 0);
    int sndbuf = 8 << 20;
    setsockopt(s, SOL_SOCKET, SO_SNDBUF, (const char *)&sndbuf, sizeof sndbuf);
    struct sockaddr_in to;
    memset(&to, 0, sizeof to);
    to.sin_family = AF_INET;
    to.sin_port = htons((unsigned short)port);
    inet_pton(AF_INET, host, &to.sin_addr);

    unsigned char *buf = malloc(65536 + 4);
    uint64_t sent = 0, pris = 0;
    uint32_t pulseOffset = 0, lastPulse = 0;
    const double tStart = now_seconds();
    double t0File = -1, loopShift = 0, lastTime = 0;
    for (int loop = 0; loop < loops; loop++) {
        FILE *f = fopen(file, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", file); return 1; }
        uint32_t currentPulse = 0xFFFFFFFFu;
        unsigned char len4[4];
        while (fread(len4, 1, 4, f) == 4) {
            const uint32_t len = get_u32(len4);
            if (len > 65536 || fread(buf, 1, len, f) != len) break;
            if (len >= 64 && get_u32(buf) == 0x52444443u) {
                const uint32_t pulse = get_u32(buf + 4);
                if (pulse != currentPulse) {                       /* first datagram of a PRI */
                    currentPulse = pulse;
                    double t;
                    memcpy(&t, buf + 32, sizeof t);
                    if (t0File < 0) t0File = t;
                    lastTime = t + loopShift;
                    if (speed > 0) sleep_until(tStart + (lastTime - t0File) / speed);
                    pris++;
                }
                lastPulse = pulse + pulseOffset;
                put_u32(buf + 4, lastPulse);                        /* continue the numbers */
                {
                    const double tShift = lastTime;                /* continue the time */
                    memcpy(buf + 32, &tShift, sizeof tShift);
                }
            }
            sendto(s, (const char *)buf, (int)len, 0, (struct sockaddr *)&to, sizeof to);
            sent++;
        }
        fclose(f);
        pulseOffset = lastPulse;                                    /* next loop: pulse numbers go on */
        loopShift = lastTime - t0File + 1e-3;
    }
    const double dt = now_seconds() - tStart;
    printf("sent %llu datagrams, %llu PRIs in %.2f s (%.0f PRIs/s)\n", (unsigned long long)sent,
           (unsigned long long)pris, dt, pris / (dt > 0 ? dt : 1));
    free(buf);
    return 0;
}
