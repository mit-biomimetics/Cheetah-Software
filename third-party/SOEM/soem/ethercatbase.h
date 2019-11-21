/*
 * Simple Open EtherCAT Master Library
 *
 * File    : ethercatbase.h
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
 * Headerfile for ethercatbase.c
 */

#ifndef _ethercatbase_
#define _ethercatbase_

#ifdef __cplusplus
extern "C"
{
#endif

int ecx_setupdatagram(ecx_portt *port, void *frame, uint8 com, uint8 idx, uint16 ADP, uint16 ADO, uint16 length, void *data);
int ecx_adddatagram(ecx_portt *port, void *frame, uint8 com, uint8 idx, boolean more, uint16 ADP, uint16 ADO, uint16 length, void *data);
int ecx_BWR(ecx_portt *port, uint16 ADP,uint16 ADO,uint16 length,void *data,int timeout);
int ecx_BRD(ecx_portt *port, uint16 ADP,uint16 ADO,uint16 length,void *data,int timeout);
int ecx_APRD(ecx_portt *port, uint16 ADP, uint16 ADO, uint16 length, void *data, int timeout);
int ecx_ARMW(ecx_portt *port, uint16 ADP, uint16 ADO, uint16 length, void *data, int timeout);
int ecx_FRMW(ecx_portt *port, uint16 ADP, uint16 ADO, uint16 length, void *data, int timeout);
uint16 ecx_APRDw(ecx_portt *port, uint16 ADP, uint16 ADO, int timeout);
int ecx_FPRD(ecx_portt *port, uint16 ADP, uint16 ADO, uint16 length, void *data, int timeout);
uint16 ecx_FPRDw(ecx_portt *port, uint16 ADP, uint16 ADO, int timeout);
int ecx_APWRw(ecx_portt *port, uint16 ADP, uint16 ADO, uint16 data, int timeout);
int ecx_APWR(ecx_portt *port, uint16 ADP, uint16 ADO, uint16 length, void *data, int timeout);
int ecx_FPWRw(ecx_portt *port, uint16 ADP, uint16 ADO, uint16 data, int timeout);
int ecx_FPWR(ecx_portt *port, uint16 ADP, uint16 ADO, uint16 length, void *data, int timeout);
int ecx_LRW(ecx_portt *port, uint32 LogAdr, uint16 length, void *data, int timeout);
int ecx_LRD(ecx_portt *port, uint32 LogAdr, uint16 length, void *data, int timeout);
int ecx_LWR(ecx_portt *port, uint32 LogAdr, uint16 length, void *data, int timeout);
int ecx_LRWDC(ecx_portt *port, uint32 LogAdr, uint16 length, void *data, uint16 DCrs, int64 *DCtime, int timeout);

#ifdef EC_VER1
int ec_setupdatagram(void *frame, uint8 com, uint8 idx, uint16 ADP, uint16 ADO, uint16 length, void *data);
int ec_adddatagram(void *frame, uint8 com, uint8 idx, boolean more, uint16 ADP, uint16 ADO, uint16 length, void *data);
int ec_BWR(uint16 ADP,uint16 ADO,uint16 length,void *data,int timeout);
int ec_BRD(uint16 ADP,uint16 ADO,uint16 length,void *data,int timeout);
int ec_APRD(uint16 ADP, uint16 ADO, uint16 length, void *data, int timeout);
int ec_ARMW(uint16 ADP, uint16 ADO, uint16 length, void *data, int timeout);
int ec_FRMW(uint16 ADP, uint16 ADO, uint16 length, void *data, int timeout);
uint16 ec_APRDw(uint16 ADP, uint16 ADO, int timeout);
int ec_FPRD(uint16 ADP, uint16 ADO, uint16 length, void *data, int timeout);
uint16 ec_FPRDw(uint16 ADP, uint16 ADO, int timeout);
int ec_APWRw(uint16 ADP, uint16 ADO, uint16 data, int timeout);
int ec_APWR(uint16 ADP, uint16 ADO, uint16 length, void *data, int timeout);
int ec_FPWRw(uint16 ADP, uint16 ADO, uint16 data, int timeout);
int ec_FPWR(uint16 ADP, uint16 ADO, uint16 length, void *data, int timeout);
int ec_LRW(uint32 LogAdr, uint16 length, void *data, int timeout);
int ec_LRD(uint32 LogAdr, uint16 length, void *data, int timeout);
int ec_LWR(uint32 LogAdr, uint16 length, void *data, int timeout);
int ec_LRWDC(uint32 LogAdr, uint16 length, void *data, uint16 DCrs, int64 *DCtime, int timeout);
#endif

#ifdef __cplusplus
}
#endif

#endif
