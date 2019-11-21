/*
 * Simple Open EtherCAT Master Library
 *
 * File    : nicdrv.h
 * Version : 1.3.1
 * Date    : 11-03-2015
 * Copyright (C) 2005-2015 Speciaal Machinefabriek Ketels v.o.f.
 * Copyright (C) 2005-2015 Arthur Ketels
 * Copyright (C) 2008-2009 TU/e Technische Universiteit Eindhoven
 * Copyright (C) 2014-2015 rt-labs AB , Sweden
 *
 * SOEM is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 2 as published by the Free
 * Software Foundation.
 *
 * SOEM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * As a special exception, if other files instantiate templates or use macros
 * or inline functions from this file, or you compile this file and link it
 * with other works to produce a work based on this file, this file does not
 * by itself cause the resulting work to be covered by the GNU General Public
 * License. However the source code for this file must still be made available
 * in accordance with section (3) of the GNU General Public License.
 *
 * This exception does not invalidate any other reasons why a work based on
 * this file might be covered by the GNU General Public License.
 *
 * The EtherCAT Technology, the trade name and logo “EtherCAT” are the intellectual
 * property of, and protected by Beckhoff Automation GmbH. You can use SOEM for
 * the sole purpose of creating, using and/or selling or otherwise distributing
 * an EtherCAT network master provided that an EtherCAT Master License is obtained
 * from Beckhoff Automation GmbH.
 *
 * In case you did not receive a copy of the EtherCAT Master License along with
 * SOEM write to Beckhoff Automation GmbH, Eiserstraße 5, D-33415 Verl, Germany
 * (www.beckhoff.com).
 */

/** \file
 * \brief
 * Headerfile for nicdrv.c
 */

#ifndef _nicdrvh_
#define _nicdrvh_

#ifdef __cplusplus
extern "C"
{
#endif

#include <pthread.h>

/** pointer structure to Tx and Rx stacks */
typedef struct
{
   /** socket connection used */
   int         *sock;
   /** tx buffer */
   ec_bufT     (*txbuf)[EC_MAXBUF];
   /** tx buffer lengths */
   int         (*txbuflength)[EC_MAXBUF];
   /** temporary receive buffer */
   ec_bufT     *tempbuf;
   /** rx buffers */
   ec_bufT     (*rxbuf)[EC_MAXBUF];
   /** rx buffer status fields */
   int         (*rxbufstat)[EC_MAXBUF];
   /** received MAC source address (middle word) */
   int         (*rxsa)[EC_MAXBUF];
} ec_stackT;

/** pointer structure to buffers for redundant port */
typedef struct
{
   ec_stackT   stack;
   int         sockhandle;
   /** rx buffers */
   ec_bufT rxbuf[EC_MAXBUF];
   /** rx buffer status */
   int rxbufstat[EC_MAXBUF];
   /** rx MAC source address */
   int rxsa[EC_MAXBUF];
   /** temporary rx buffer */
   ec_bufT tempinbuf;
} ecx_redportt;

/** pointer structure to buffers, vars and mutexes for port instantiation */
typedef struct
{
   ec_stackT   stack;
   int         sockhandle;
   /** rx buffers */
   ec_bufT rxbuf[EC_MAXBUF];
   /** rx buffer status */
   int rxbufstat[EC_MAXBUF];
   /** rx MAC source address */
   int rxsa[EC_MAXBUF];
   /** temporary rx buffer */
   ec_bufT tempinbuf;
   /** temporary rx buffer status */
   int tempinbufs;
   /** transmit buffers */
   ec_bufT txbuf[EC_MAXBUF];
   /** transmit buffer lenghts */
   int txbuflength[EC_MAXBUF];
   /** temporary tx buffer */
   ec_bufT txbuf2;
   /** temporary tx buffer length */
   int txbuflength2;
   /** last used frame index */
   int lastidx;
   /** current redundancy state */
   int redstate;
   /** pointer to redundancy port and buffers */
   ecx_redportt *redport;
   pthread_mutex_t getindex_mutex;
   pthread_mutex_t tx_mutex;
   pthread_mutex_t rx_mutex;
} ecx_portt;

extern const uint16 priMAC[3];
extern const uint16 secMAC[3];

#ifdef EC_VER1
extern ecx_portt     ecx_port;
extern ecx_redportt  ecx_redport;

int ec_setupnic(const char * ifname, int secondary);
int ec_closenic(void);
void ec_setbufstat(int idx, int bufstat);
int ec_getindex(void);
int ec_outframe(int idx, int sock);
int ec_outframe_red(int idx);
int ec_waitinframe(int idx, int timeout);
int ec_srconfirm(int idx,int timeout);
#endif

void ec_setupheader(void *p);
int ecx_setupnic(ecx_portt *port, const char * ifname, int secondary);
int ecx_closenic(ecx_portt *port);
void ecx_setbufstat(ecx_portt *port, int idx, int bufstat);
int ecx_getindex(ecx_portt *port);
int ecx_outframe(ecx_portt *port, int idx, int sock);
int ecx_outframe_red(ecx_portt *port, int idx);
int ecx_waitinframe(ecx_portt *port, int idx, int timeout);
int ecx_srconfirm(ecx_portt *port, int idx,int timeout);

#ifdef __cplusplus
}
#endif

#endif
