/*
 * Simple Open EtherCAT Master Library
 *
 * File    : ethercatdc.c
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
 * Distributed Clock EtherCAT functions.
 *
 */
#include "oshw.h"
#include "osal.h"
#include "ethercattype.h"
#include "ethercatbase.h"
#include "ethercatmain.h"
#include "ethercatdc.h"

#define PORTM0 0x01
#define PORTM1 0x02
#define PORTM2 0x04
#define PORTM3 0x08

/** 1st sync pulse delay in ns here 100ms */
#define SyncDelay       ((int32)100000000)

/**
 * Set DC of slave to fire sync0 at CyclTime interval with CyclShift offset.
 *
 * @param[in]  context        = context struct
 * @param [in] slave            Slave number.
 * @param [in] act              TRUE = active, FALSE = deactivated
 * @param [in] CyclTime         Cycltime in ns.
 * @param [in] CyclShift        CyclShift in ns.
 */
void ecx_dcsync0(ecx_contextt *context, uint16 slave, boolean act, uint32 CyclTime, int32 CyclShift)
{
   uint8 h, RA;
   uint16 slaveh;
   int64 t, t1;
   int32 tc;

   slaveh = context->slavelist[slave].configadr;
   RA = 0;

   /* stop cyclic operation, ready for next trigger */
   (void)ecx_FPWR(context->port, slaveh, ECT_REG_DCSYNCACT, sizeof(RA), &RA, EC_TIMEOUTRET);
   if (act)
   {
       RA = 1 + 2;    /* act cyclic operation and sync0, sync1 deactivated */
   }
   h = 0;
   (void)ecx_FPWR(context->port, slaveh, ECT_REG_DCCUC, sizeof(h), &h, EC_TIMEOUTRET); /* write access to ethercat */
   t1 = 0;
   (void)ecx_FPRD(context->port, slaveh, ECT_REG_DCSYSTIME, sizeof(t1), &t1, EC_TIMEOUTRET); /* read local time of slave */
   t1 = etohll(t1);

   /* Calculate first trigger time, always a whole multiple of CyclTime rounded up
   plus the shifttime (can be negative)
   This insures best sychronisation between slaves, slaves with the same CyclTime
   will sync at the same moment (you can use CyclShift to shift the sync) */
   if (CyclTime > 0)
   {
       t = ((t1 + SyncDelay) / CyclTime) * CyclTime + CyclTime + CyclShift;
   }
   else
   {
      t = t1 + SyncDelay + CyclShift;
      /* first trigger at T1 + CyclTime + SyncDelay + CyclShift in ns */
   }
   t = htoell(t);
   (void)ecx_FPWR(context->port, slaveh, ECT_REG_DCSTART0, sizeof(t), &t, EC_TIMEOUTRET); /* SYNC0 start time */
   tc = htoel(CyclTime);
   (void)ecx_FPWR(context->port, slaveh, ECT_REG_DCCYCLE0, sizeof(tc), &tc, EC_TIMEOUTRET); /* SYNC0 cycle time */
   (void)ecx_FPWR(context->port, slaveh, ECT_REG_DCSYNCACT, sizeof(RA), &RA, EC_TIMEOUTRET); /* activate cyclic operation */

    // update ec_slave state
    context->slavelist[slave].DCactive = (uint8)act;
    context->slavelist[slave].DCshift = CyclShift;
    context->slavelist[slave].DCcycle = CyclTime;
}

/**
 * Set DC of slave to fire sync0 and sync1 at CyclTime interval with CyclShift offset.
 *
 * @param[in]  context        = context struct
 * @param [in] slave            Slave number.
 * @param [in] act              TRUE = active, FALSE = deactivated
 * @param [in] CyclTime0        Cycltime SYNC0 in ns.
 * @param [in] CyclTime1        Cycltime SYNC1 in ns. This time is a delta time in relation to
                                the SYNC0 fire. If CylcTime1 = 0 then SYNC1 fires a the same time
                                as SYNC0.
 * @param [in] CyclShift        CyclShift in ns.
 */
void ecx_dcsync01(ecx_contextt *context, uint16 slave, boolean act, uint32 CyclTime0, uint32 CyclTime1, int32 CyclShift)
{
   uint8 h, RA;
   uint16 slaveh;
   int64 t, t1;
   int32 tc;
   uint32 TrueCyclTime;

   /* Sync1 can be used as a multiple of Sync0, use true cycle time */
   TrueCyclTime = ((CyclTime1 / CyclTime0) + 1) * CyclTime0;

   slaveh = context->slavelist[slave].configadr;
   RA = 0;

   /* stop cyclic operation, ready for next trigger */
   (void)ecx_FPWR(context->port, slaveh, ECT_REG_DCSYNCACT, sizeof(RA), &RA, EC_TIMEOUTRET);
   if (act)
   {
      RA = 1 + 2 + 4;    /* act cyclic operation and sync0 + sync1 */
   }
   h = 0;
   (void)ecx_FPWR(context->port, slaveh, ECT_REG_DCCUC, sizeof(h), &h, EC_TIMEOUTRET); /* write access to ethercat */
   t1 = 0;
   (void)ecx_FPRD(context->port, slaveh, ECT_REG_DCSYSTIME, sizeof(t1), &t1, EC_TIMEOUTRET); /* read local time of slave */
   t1 = etohll(t1);

   /* Calculate first trigger time, always a whole multiple of TrueCyclTime rounded up
   plus the shifttime (can be negative)
   This insures best sychronisation between slaves, slaves with the same CyclTime
   will sync at the same moment (you can use CyclShift to shift the sync) */
   if (CyclTime0 > 0)
   {
      t = ((t1 + SyncDelay) / TrueCyclTime) * TrueCyclTime + TrueCyclTime + CyclShift;
   }
   else
   {
      t = t1 + SyncDelay + CyclShift;
      /* first trigger at T1 + CyclTime + SyncDelay + CyclShift in ns */
   }
   t = htoell(t);
   (void)ecx_FPWR(context->port, slaveh, ECT_REG_DCSTART0, sizeof(t), &t, EC_TIMEOUTRET); /* SYNC0 start time */
   tc = htoel(CyclTime0);
   (void)ecx_FPWR(context->port, slaveh, ECT_REG_DCCYCLE0, sizeof(tc), &tc, EC_TIMEOUTRET); /* SYNC0 cycle time */
   tc = htoel(CyclTime1);
   (void)ecx_FPWR(context->port, slaveh, ECT_REG_DCCYCLE1, sizeof(tc), &tc, EC_TIMEOUTRET); /* SYNC1 cycle time */
   (void)ecx_FPWR(context->port, slaveh, ECT_REG_DCSYNCACT, sizeof(RA), &RA, EC_TIMEOUTRET); /* activate cyclic operation */

    // update ec_slave state
    context->slavelist[slave].DCactive = (uint8)act;
    context->slavelist[slave].DCshift = CyclShift;
    context->slavelist[slave].DCcycle = CyclTime0;
}

/* latched port time of slave */
static int32 ecx_porttime(ecx_contextt *context, uint16 slave, uint8 port)
{
   int32 ts;
   switch (port)
   {
      case 0:
         ts = context->slavelist[slave].DCrtA;
         break;
      case 1:
         ts = context->slavelist[slave].DCrtB;
         break;
      case 2:
         ts = context->slavelist[slave].DCrtC;
         break;
      case 3:
         ts = context->slavelist[slave].DCrtD;
         break;
      default:
         ts = 0;
         break;
   }
   return ts;
}

/* calculate previous active port of a slave */
static uint8 ecx_prevport(ecx_contextt *context, uint16 slave, uint8 port)
{
   uint8 pport = port;
   uint8 aport = context->slavelist[slave].activeports;
   switch(port)
   {
      case 0:
         if(aport & PORTM2)
            pport = 2;
         else if (aport & PORTM1)
            pport = 1;
         else if (aport & PORTM3)
            pport = 3;
         break;
      case 1:
         if(aport & PORTM3)
            pport = 3;
         else if (aport & PORTM0)
            pport = 0;
         else if (aport & PORTM2)
            pport = 2;
         break;
      case 2:
         if(aport & PORTM1)
            pport = 1;
         else if (aport & PORTM3)
            pport = 3;
         else if (aport & PORTM0)
            pport = 0;
         break;
      case 3:
         if(aport & PORTM0)
            pport = 0;
         else if (aport & PORTM2)
            pport = 2;
         else if (aport & PORTM1)
            pport = 1;
         break;
   }
   return pport;
}

/* search unconsumed ports in parent, consume and return first open port */
static uint8 ecx_parentport(ecx_contextt *context, uint16 parent)
{
   uint8 parentport = 0;
   uint8 b;
   /* search order is important, here 3 - 1 - 2 - 0 */
   b = context->slavelist[parent].consumedports;
   if (b & PORTM3)
   {
      parentport = 3;
      b &= (uint8)~PORTM3;
   }
   else if (b & PORTM1)
   {
      parentport = 1;
      b &= (uint8)~PORTM1;
   }
   else if (b & PORTM2)
   {
      parentport = 2;
      b &= (uint8)~PORTM2;
   }
   else if (b & PORTM0)
   {
      parentport = 0;
      b &= (uint8)~PORTM0;
   }
   context->slavelist[parent].consumedports = b;
   return parentport;
}

/**
 * Locate DC slaves, measure propagation delays.
 *
 * @param[in]  context        = context struct
 * @return boolean if slaves are found with DC
 */
boolean ecx_configdc(ecx_contextt *context)
{
   uint16 i, slaveh, parent, child;
   uint16 parenthold = 0;
   uint16 prevDCslave = 0;
   int32 ht, dt1, dt2, dt3;
   int64 hrt;
   uint8 entryport;
   int8 nlist;
   int8 plist[4];
   int32 tlist[4];
   ec_timet mastertime;
   uint64 mastertime64;

   context->slavelist[0].hasdc = FALSE;
   context->grouplist[0].hasdc = FALSE;
   ht = 0;

   ecx_BWR(context->port, 0, ECT_REG_DCTIME0, sizeof(ht), &ht, EC_TIMEOUTRET);  /* latch DCrecvTimeA of all slaves */
   mastertime = osal_current_time();
   mastertime.sec -= 946684800UL;  /* EtherCAT uses 2000-01-01 as epoch start instead of 1970-01-01 */
   mastertime64 = (((uint64)mastertime.sec * 1000000) + (uint64)mastertime.usec) * 1000;
   for (i = 1; i <= *(context->slavecount); i++)
   {
      context->slavelist[i].consumedports = context->slavelist[i].activeports;
      if (context->slavelist[i].hasdc)
      {
         if (!context->slavelist[0].hasdc)
         {
            context->slavelist[0].hasdc = TRUE;
            context->slavelist[0].DCnext = i;
            context->slavelist[i].DCprevious = 0;
            context->grouplist[0].hasdc = TRUE;
            context->grouplist[0].DCnext = i;
         }
         else
         {
            context->slavelist[prevDCslave].DCnext = i;
            context->slavelist[i].DCprevious = prevDCslave;
         }
         /* this branch has DC slave so remove parenthold */
         parenthold = 0;
         prevDCslave = i;
         slaveh = context->slavelist[i].configadr;
         (void)ecx_FPRD(context->port, slaveh, ECT_REG_DCTIME0, sizeof(ht), &ht, EC_TIMEOUTRET);
         context->slavelist[i].DCrtA = etohl(ht);
         /* 64bit latched DCrecvTimeA of each specific slave */
         (void)ecx_FPRD(context->port, slaveh, ECT_REG_DCSOF, sizeof(hrt), &hrt, EC_TIMEOUTRET);
         /* use it as offset in order to set local time around 0 + mastertime */
         hrt = htoell(-etohll(hrt) + mastertime64);
         /* save it in the offset register */
         (void)ecx_FPWR(context->port, slaveh, ECT_REG_DCSYSOFFSET, sizeof(hrt), &hrt, EC_TIMEOUTRET);
         (void)ecx_FPRD(context->port, slaveh, ECT_REG_DCTIME1, sizeof(ht), &ht, EC_TIMEOUTRET);
         context->slavelist[i].DCrtB = etohl(ht);
         (void)ecx_FPRD(context->port, slaveh, ECT_REG_DCTIME2, sizeof(ht), &ht, EC_TIMEOUTRET);
         context->slavelist[i].DCrtC = etohl(ht);
         (void)ecx_FPRD(context->port, slaveh, ECT_REG_DCTIME3, sizeof(ht), &ht, EC_TIMEOUTRET);
         context->slavelist[i].DCrtD = etohl(ht);

         /* make list of active ports and their time stamps */
         nlist = 0;
         if (context->slavelist[i].activeports & PORTM0)
         {
            plist[nlist] = 0;
            tlist[nlist] = context->slavelist[i].DCrtA;
            nlist++;
         }
         if (context->slavelist[i].activeports & PORTM3)
         {
            plist[nlist] = 3;
            tlist[nlist] = context->slavelist[i].DCrtD;
            nlist++;
         }
         if (context->slavelist[i].activeports & PORTM1)
         {
            plist[nlist] = 1;
            tlist[nlist] = context->slavelist[i].DCrtB;
            nlist++;
         }
         if (context->slavelist[i].activeports & PORTM2)
         {
            plist[nlist] = 2;
            tlist[nlist] = context->slavelist[i].DCrtC;
            nlist++;
         }
         /* entryport is port with the lowest timestamp */
         entryport = 0;
         if((nlist > 1) && (tlist[1] < tlist[entryport]))
         {
            entryport = 1;
         }
         if((nlist > 2) && (tlist[2] < tlist[entryport]))
         {
            entryport = 2;
         }
         if((nlist > 3) && (tlist[3] < tlist[entryport]))
         {
            entryport = 3;
         }
         entryport = plist[entryport];
         context->slavelist[i].entryport = entryport;
         /* consume entryport from activeports */
         context->slavelist[i].consumedports &= (uint8)~(1 << entryport);

         /* finding DC parent of current */
         parent = i;
         do
         {
            child = parent;
            parent = context->slavelist[parent].parent;
         }
         while (!((parent == 0) || (context->slavelist[parent].hasdc)));
         /* only calculate propagation delay if slave is not the first */
         if (parent > 0)
         {
            /* find port on parent this slave is connected to */
            context->slavelist[i].parentport = ecx_parentport(context, parent);
            if (context->slavelist[parent].topology == 1)
            {
               context->slavelist[i].parentport = context->slavelist[parent].entryport;
            }

            dt1 = 0;
            dt2 = 0;
            /* delta time of (parentport - 1) - parentport */
            /* note: order of ports is 0 - 3 - 1 -2 */
            /* non active ports are skipped */
            dt3 = ecx_porttime(context, parent, context->slavelist[i].parentport) -
                  ecx_porttime(context, parent,
                    ecx_prevport(context, parent, context->slavelist[i].parentport));
            /* current slave has children */
            /* those children's delays need to be subtracted */
            if (context->slavelist[i].topology > 1)
            {
               dt1 = ecx_porttime(context, i,
                        ecx_prevport(context, i, context->slavelist[i].entryport)) -
                     ecx_porttime(context, i, context->slavelist[i].entryport);
            }
            /* we are only interested in positive difference */
            if (dt1 > dt3) dt1 = -dt1;
            /* current slave is not the first child of parent */
            /* previous child's delays need to be added */
            if ((child - parent) > 1)
            {
               dt2 = ecx_porttime(context, parent,
                        ecx_prevport(context, parent, context->slavelist[i].parentport)) -
                     ecx_porttime(context, parent, context->slavelist[parent].entryport);
            }
            if (dt2 < 0) dt2 = -dt2;

            /* calculate current slave delay from delta times */
            /* assumption : forward delay equals return delay */
            context->slavelist[i].pdelay = ((dt3 - dt1) / 2) + dt2 +
               context->slavelist[parent].pdelay;
            ht = htoel(context->slavelist[i].pdelay);
            /* write propagation delay*/
            (void)ecx_FPWR(context->port, slaveh, ECT_REG_DCSYSDELAY, sizeof(ht), &ht, EC_TIMEOUTRET);
         }
      }
      else
      {
         context->slavelist[i].DCrtA = 0;
         context->slavelist[i].DCrtB = 0;
         context->slavelist[i].DCrtC = 0;
         context->slavelist[i].DCrtD = 0;
         parent = context->slavelist[i].parent;
         /* if non DC slave found on first position on branch hold root parent */
         if ( (parent > 0) && (context->slavelist[parent].topology > 2))
            parenthold = parent;
         /* if branch has no DC slaves consume port on root parent */
         if ( parenthold && (context->slavelist[i].topology == 1))
         {
            ecx_parentport(context, parenthold);
            parenthold = 0;
         }
      }
   }

   return context->slavelist[0].hasdc;
}

#ifdef EC_VER1
void ec_dcsync0(uint16 slave, boolean act, uint32 CyclTime, int32 CyclShift)
{
   ecx_dcsync0(&ecx_context, slave, act, CyclTime, CyclShift);
}

void ec_dcsync01(uint16 slave, boolean act, uint32 CyclTime0, uint32 CyclTime1, int32 CyclShift)
{
   ecx_dcsync01(&ecx_context, slave, act, CyclTime0, CyclTime1, CyclShift);
}

boolean ec_configdc(void)
{
   return ecx_configdc(&ecx_context);
}
#endif
