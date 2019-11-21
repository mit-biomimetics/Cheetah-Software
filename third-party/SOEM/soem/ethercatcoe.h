/*
 * Simple Open EtherCAT Master Library
 *
 * File    : ethercatcoe.h
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
 * Headerfile for ethercatcoe.c
 */

#ifndef _ethercatcoe_
#define _ethercatcoe_

#ifdef __cplusplus
extern "C"
{
#endif

/** max entries in Object Description list */
#define EC_MAXODLIST   1024

/** max entries in Object Entry list */
#define EC_MAXOELIST   256

/* Storage for object description list */
typedef struct
{
   /** slave number */
   uint16  Slave;
   /** number of entries in list */
   uint16  Entries;
   /** array of indexes */
   uint16  Index[EC_MAXODLIST];
   /** array of datatypes, see EtherCAT specification */
   uint16  DataType[EC_MAXODLIST];
   /** array of object codes, see EtherCAT specification */
   uint8   ObjectCode[EC_MAXODLIST];
   /** number of subindexes for each index */
   uint8   MaxSub[EC_MAXODLIST];
   /** textual description of each index */
   char    Name[EC_MAXODLIST][EC_MAXNAME+1];
} ec_ODlistt;

/* storage for object list entry information */
typedef struct
{
   /** number of entries in list */
   uint16 Entries;
   /** array of value infos, see EtherCAT specification */
   uint8  ValueInfo[EC_MAXOELIST];
   /** array of value infos, see EtherCAT specification */
   uint16 DataType[EC_MAXOELIST];
   /** array of bit lengths, see EtherCAT specification */
   uint16 BitLength[EC_MAXOELIST];
   /** array of object access bits, see EtherCAT specification */
   uint16 ObjAccess[EC_MAXOELIST];
   /** textual description of each index */
   char   Name[EC_MAXOELIST][EC_MAXNAME+1];
} ec_OElistt;

#ifdef EC_VER1
void ec_SDOerror(uint16 Slave, uint16 Index, uint8 SubIdx, int32 AbortCode);
int ec_SDOread(uint16 slave, uint16 index, uint8 subindex,
                      boolean CA, int *psize, void *p, int timeout);
int ec_SDOwrite(uint16 Slave, uint16 Index, uint8 SubIndex,
    boolean CA, int psize, void *p, int Timeout);
int ec_RxPDO(uint16 Slave, uint16 RxPDOnumber , int psize, void *p);
int ec_TxPDO(uint16 slave, uint16 TxPDOnumber , int *psize, void *p, int timeout);
int ec_readPDOmap(uint16 Slave, int *Osize, int *Isize);
int ec_readPDOmapCA(uint16 Slave, int Thread_n, int *Osize, int *Isize);
int ec_readODlist(uint16 Slave, ec_ODlistt *pODlist);
int ec_readODdescription(uint16 Item, ec_ODlistt *pODlist);
int ec_readOEsingle(uint16 Item, uint8 SubI, ec_ODlistt *pODlist, ec_OElistt *pOElist);
int ec_readOE(uint16 Item, ec_ODlistt *pODlist, ec_OElistt *pOElist);
#endif

void ecx_SDOerror(ecx_contextt *context, uint16 Slave, uint16 Index, uint8 SubIdx, int32 AbortCode);
int ecx_SDOread(ecx_contextt *context, uint16 slave, uint16 index, uint8 subindex,
                      boolean CA, int *psize, void *p, int timeout);
int ecx_SDOwrite(ecx_contextt *context, uint16 Slave, uint16 Index, uint8 SubIndex,
    boolean CA, int psize, void *p, int Timeout);
int ecx_RxPDO(ecx_contextt *context, uint16 Slave, uint16 RxPDOnumber , int psize, void *p);
int ecx_TxPDO(ecx_contextt *context, uint16 slave, uint16 TxPDOnumber , int *psize, void *p, int timeout);
int ecx_readPDOmap(ecx_contextt *context, uint16 Slave, int *Osize, int *Isize);
int ecx_readPDOmapCA(ecx_contextt *context, uint16 Slave, int Thread_n, int *Osize, int *Isize);
int ecx_readODlist(ecx_contextt *context, uint16 Slave, ec_ODlistt *pODlist);
int ecx_readODdescription(ecx_contextt *context, uint16 Item, ec_ODlistt *pODlist);
int ecx_readOEsingle(ecx_contextt *context, uint16 Item, uint8 SubI, ec_ODlistt *pODlist, ec_OElistt *pOElist);
int ecx_readOE(ecx_contextt *context, uint16 Item, ec_ODlistt *pODlist, ec_OElistt *pOElist);

#ifdef __cplusplus
}
#endif

#endif
