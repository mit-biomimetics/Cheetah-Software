/*!
 * @file rt_serial.cpp
 * @brief Serial port
 */

#ifdef linux

#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define termios asmtermios

#include <asm/termios.h>

#undef termios

#include <termios.h>
#include <math.h>
#include <pthread.h>
#include <stropts.h>
#include <endian.h>
#include <stdint.h>

#include "rt/rt_serial.h"

void init_serial_for_sbus(int fd, int baud) {
  printf("\t[RT SERIAL] Configuring serial device...\n");
  struct termios2 tty;

  ioctl(fd, TCGETS2, &tty);
  tty.c_cflag &= ~CBAUD;
  tty.c_cflag |= BOTHER;
  tty.c_ispeed = baud;
  tty.c_ospeed = baud;


tty.c_iflag &= ~(IGNBRK | BRKINT | PARMRK | ISTRIP
                | INLCR | IGNCR | ICRNL | IXON);
tty.c_oflag &= ~OPOST;
tty.c_lflag &= ~(ECHO | ECHONL | ICANON | ISIG | IEXTEN);
tty.c_cflag &= ~(CSIZE | PARENB);
tty.c_cflag |= PARENB;
tty.c_cflag &= ~PARODD;
tty.c_cflag |=
tty.c_cflag |= CS8;

  // tty.c_cflag = (tty.c_cflag & ~CSIZE) | CS8;  // 8-bit chars
  // disable IGNBRK for mismatched speed tests; otherwise receive break
  // as \000 chars
  // tty.c_iflag &= ~IGNBRK;  // disable break processing
  // tty.c_lflag = 0;         // no signaling chars, no echo,
  // // no canonical processing
  // tty.c_oflag = 0;      // no remapping, no delays
  tty.c_cc[VMIN] = 1;   // read doesn't block
  tty.c_cc[VTIME] = 1;  // 0.5 seconds read timeout

  // tty.c_iflag &= ~(IXON | IXOFF | IXANY);  // shut off xon/xoff ctrl

  // tty.c_cflag |= (CLOCAL | CREAD);  // ignore modem controls,
  // // enable reading
  // // tty.c_cflag &= ~(PARENB | PARODD);      // shut off parity
  // tty.c_cflag |= PARENB;
  // tty.c_cflag &= ~CSTOPB;
  // tty.c_cflag &= ~CRTSCTS;
  // cfmakeraw(&tty);

  ioctl(fd, TCSETS2, &tty);
}

/**
 * @brief Configure serial port
 *
 * @param fd Serial port FD
 * @param speed Baud rate
 * @param parity Parity
 * @param port Port number
 */
int set_interface_attribs_custom_baud(int fd, int speed, int parity, int port) {
  (void)parity;
  (void)port;

  printf("\t[RT SERIAL] Configuring serial device...\n");
  struct termios2 tty;

  ioctl(fd, TCGETS2, &tty);
  tty.c_cflag &= ~CBAUD;
  tty.c_cflag |= BOTHER;
  tty.c_ispeed = speed;
  tty.c_ospeed = speed;

  tty.c_cflag = (tty.c_cflag & ~CSIZE) | CS8;  // 8-bit chars
  // disable IGNBRK for mismatched speed tests; otherwise receive break
  // as \000 chars
  tty.c_iflag &= ~IGNBRK;  // disable break processing
  tty.c_lflag = 0;         // no signaling chars, no echo,
  // no canonical processing
  tty.c_oflag = 0;      // no remapping, no delays
  tty.c_cc[VMIN] = 0;   // read doesn't block
  tty.c_cc[VTIME] = 1;  // 0.5 seconds read timeout

  tty.c_iflag &= ~(IXON | IXOFF | IXANY);  // shut off xon/xoff ctrl

  tty.c_cflag |= (CLOCAL | CREAD);  // ignore modem controls,
  // enable reading
  // tty.c_cflag &= ~(PARENB | PARODD);      // shut off parity
  tty.c_cflag |= PARENB;
  tty.c_cflag &= ~CSTOPB;
  tty.c_cflag &= ~CRTSCTS;
  // cfmakeraw(&tty);

  ioctl(fd, TCSETS2, &tty);
  return 0;
}

#endif
