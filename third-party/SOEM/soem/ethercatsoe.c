/*
 * Simple Open EtherCAT Master Library
 *
 * File    : ethercatsoe.c
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
 * Servo over EtherCAT (SoE) Module.
 */

#include <stdio.h>
#include <string.h>
#include "osal.h"
#include "oshw.h"
#include "ethercattype.h"
#include "ethercatbase.h"
#include "ethercatmain.h"
#include "ethercatsoe.h"

#define EC_SOE_MAX_DRIVES 8

/** SoE (Servo over EtherCAT) mailbox structure */
PACKED_BEGIN
typedef struct PACKED
{
   ec_mbxheadert MbxHeader;
   uint8         opCode         :3;
   uint8         incomplete     :1;
   uint8         error          :1;
   uint8         driveNo        :3;
   uint8         elementflags;
   union
   {
      uint16     idn;
      uint16     fragmentsleft;
   };
} ec_SoEt;
PACKED_END

/** Report SoE error.
 *
 * @param[in]  context        = context struct
 * @param[in]  Slave      = Slave number
 * @param[in]  idn        = IDN that generated error
 * @param[in]  Error      = Error code, see EtherCAT documentation for list
 */
void ecx_SoEerror(ecx_contextt *context, uint16 Slave, uint16 idn, uint16 Error)
{
   ec_errort Ec;

   memset(&Ec, 0, sizeof(Ec));
   Ec.Time = osal_current_time();
   Ec.Slave = Slave;
   Ec.Index = idn;
   Ec.SubIdx = 0;
   *(context->ecaterror) = TRUE;
   Ec.Etype = EC_ERR_TYPE_SOE_ERROR;
   Ec.ErrorCode = Error;
   ecx_pusherror(context, &Ec);
}

/** SoE read, blocking.
 *
 * The IDN object of the selected slave and DriveNo is read. If a response
 * is larger than the mailbox size then the response is segmented. The function
 * will combine all segments and copy them to the parameter buffer.
 *
 * @param[in]  context        = context struct
 * @param[in]  slave         = Slave number
 * @param[in]  driveNo       = Drive number in slave
 * @param[in]  elementflags  = Flags to select what properties of IDN are to be transfered.
 * @param[in]  idn           = IDN.
 * @param[in,out] psize      = Size in bytes of parameter buffer, returns bytes read from SoE.
 * @param[out] p             = Pointer to parameter buffer
 * @param[in]  timeout       = Timeout in us, standard is EC_TIMEOUTRXM
 * @return Workcounter from last slave response
 */
int ecx_SoEread(ecx_contextt *context, uint16 slave, uint8 driveNo, uint8 elementflags, uint16 idn, int *psize, void *p, int timeout)
{
   ec_SoEt *SoEp, *aSoEp;
   uint16 totalsize, framedatasize;
   int wkc;
   uint8 *bp;
   uint8 *mp;
   uint16 *errorcode;
   ec_mbxbuft MbxIn, MbxOut;
   uint8 cnt;
   boolean NotLast;

   ec_clearmbx(&MbxIn);
   /* Empty slave out mailbox if something is in. Timeout set to 0 */
   wkc = ecx_mbxreceive(context, slave, (ec_mbxbuft *)&MbxIn, 0);
   ec_clearmbx(&MbxOut);
   aSoEp = (ec_SoEt *)&MbxIn;
   SoEp = (ec_SoEt *)&MbxOut;
   SoEp->MbxHeader.length = htoes(sizeof(ec_SoEt) - sizeof(ec_mbxheadert));
   SoEp->MbxHeader.address = htoes(0x0000);
   SoEp->MbxHeader.priority = 0x00;
   /* get new mailbox count value, used as session handle */
   cnt = ec_nextmbxcnt(context->slavelist[slave].mbx_cnt);
   context->slavelist[slave].mbx_cnt = cnt;
   SoEp->MbxHeader.mbxtype = ECT_MBXT_SOE + (cnt << 4); /* SoE */
   SoEp->opCode = ECT_SOE_READREQ;
   SoEp->incomplete = 0;
   SoEp->error = 0;
   SoEp->driveNo = driveNo;
   SoEp->elementflags = elementflags;
   SoEp->idn = htoes(idn);
   totalsize = 0;
   bp = p;
   mp = (uint8 *)&MbxIn + sizeof(ec_SoEt);
   NotLast = TRUE;
   /* send SoE request to slave */
   wkc = ecx_mbxsend(context, slave, (ec_mbxbuft *)&MbxOut, EC_TIMEOUTTXM);
   if (wkc > 0) /* succeeded to place mailbox in slave ? */
   {
      while (NotLast)
      {
         /* clean mailboxbuffer */
         ec_clearmbx(&MbxIn);
         /* read slave response */
         wkc = ecx_mbxreceive(context, slave, (ec_mbxbuft *)&MbxIn, timeout);
         if (wkc > 0) /* succeeded to read slave response ? */
         {
            /* slave response should be SoE, ReadRes */
            if (((aSoEp->MbxHeader.mbxtype & 0x0f) == ECT_MBXT_SOE) &&
                (aSoEp->opCode == ECT_SOE_READRES) &&
                (aSoEp->error == 0) &&
                (aSoEp->driveNo == driveNo) &&
                (aSoEp->elementflags == elementflags))
            {
               framedatasize = etohs(aSoEp->MbxHeader.length) - sizeof(ec_SoEt)  + sizeof(ec_mbxheadert);
               totalsize += framedatasize;
               /* Does parameter fit in parameter buffer ? */
               if (totalsize <= *psize)
               {
                  /* copy parameter data in parameter buffer */
                  memcpy(bp, mp, framedatasize);
                  /* increment buffer pointer */
                  bp += framedatasize;
               }
               else
               {
                  framedatasize -= totalsize - *psize;
                  totalsize = *psize;
                  /* copy parameter data in parameter buffer */
                  if (framedatasize > 0) memcpy(bp, mp, framedatasize);
               }

               if (!aSoEp->incomplete)
               {
                  NotLast = FALSE;
                  *psize = totalsize;
               }
            }
            /* other slave response */
            else
            {
               NotLast = FALSE;
               if (((aSoEp->MbxHeader.mbxtype & 0x0f) == ECT_MBXT_SOE) &&
                   (aSoEp->opCode == ECT_SOE_READRES) &&
                   (aSoEp->error == 1))
               {
                  mp = (uint8 *)&MbxIn + (etohs(aSoEp->MbxHeader.length) + sizeof(ec_mbxheadert) - sizeof(uint16));
                  errorcode = (uint16 *)mp;
                  ecx_SoEerror(context, slave, idn, *errorcode);
               }
               else
               {
                  ecx_packeterror(context, slave, idn, 0, 1); /* Unexpected frame returned */
               }
               wkc = 0;
            }
         }
         else
         {
            NotLast = FALSE;
            ecx_packeterror(context, slave, idn, 0, 4); /* no response */
         }
      }
   }
   return wkc;
}

/** SoE write, blocking.
 *
 * The IDN object of the selected slave and DriveNo is written. If a response
 * is larger than the mailbox size then the response is segmented.
 *
 * @param[in]  context        = context struct
 * @param[in]  slave         = Slave number
 * @param[in]  driveNo       = Drive number in slave
 * @param[in]  elementflags  = Flags to select what properties of IDN are to be transfered.
 * @param[in]  idn           = IDN.
 * @param[in]  psize         = Size in bytes of parameter buffer.
 * @param[out] p             = Pointer to parameter buffer
 * @param[in]  timeout       = Timeout in us, standard is EC_TIMEOUTRXM
 * @return Workcounter from last slave response
 */
int ecx_SoEwrite(ecx_contextt *context, uint16 slave, uint8 driveNo, uint8 elementflags, uint16 idn, int psize, void *p, int timeout)
{
   ec_SoEt *SoEp, *aSoEp;
   uint16 framedatasize, maxdata;
   int wkc;
   uint8 *mp;
   uint8 *hp;
   uint16 *errorcode;
   ec_mbxbuft MbxIn, MbxOut;
   uint8 cnt;
   boolean NotLast;

   ec_clearmbx(&MbxIn);
   /* Empty slave out mailbox if something is in. Timeout set to 0 */
   wkc = ecx_mbxreceive(context, slave, (ec_mbxbuft *)&MbxIn, 0);
   ec_clearmbx(&MbxOut);
   aSoEp = (ec_SoEt *)&MbxIn;
   SoEp = (ec_SoEt *)&MbxOut;
   SoEp->MbxHeader.address = htoes(0x0000);
   SoEp->MbxHeader.priority = 0x00;
   SoEp->opCode = ECT_SOE_WRITEREQ;
   SoEp->error = 0;
   SoEp->driveNo = driveNo;
   SoEp->elementflags = elementflags;
   hp = p;
   mp = (uint8 *)&MbxOut + sizeof(ec_SoEt);
   maxdata = context->slavelist[slave].mbx_l - sizeof(ec_SoEt);
   NotLast = TRUE;
   while (NotLast)
   {
      framedatasize = psize;
      NotLast = FALSE;
      SoEp->idn = htoes(idn);
      SoEp->incomplete = 0;
      if (framedatasize > maxdata)
      {
         framedatasize = maxdata;  /*  segmented transfer needed  */
         NotLast = TRUE;
         SoEp->incomplete = 1;
         SoEp->fragmentsleft = psize / maxdata;
      }
      SoEp->MbxHeader.length = htoes(sizeof(ec_SoEt) - sizeof(ec_mbxheadert) + framedatasize);
      /* get new mailbox counter, used for session handle */
      cnt = ec_nextmbxcnt(context->slavelist[slave].mbx_cnt);
      context->slavelist[slave].mbx_cnt = cnt;
      SoEp->MbxHeader.mbxtype = ECT_MBXT_SOE + (cnt << 4); /* SoE */
      /* copy parameter data to mailbox */
      memcpy(mp, hp, framedatasize);
      hp += framedatasize;
      psize -= framedatasize;
      /* send SoE request to slave */
      wkc = ecx_mbxsend(context, slave, (ec_mbxbuft *)&MbxOut, EC_TIMEOUTTXM);
      if (wkc > 0) /* succeeded to place mailbox in slave ? */
      {
         if (!NotLast || !ecx_mbxempty(context, slave, timeout))
         {
            /* clean mailboxbuffer */
            ec_clearmbx(&MbxIn);
            /* read slave response */
            wkc = ecx_mbxreceive(context, slave, (ec_mbxbuft *)&MbxIn, timeout);
            if (wkc > 0) /* succeeded to read slave response ? */
            {
               NotLast = FALSE;
               /* slave response should be SoE, WriteRes */
               if (((aSoEp->MbxHeader.mbxtype & 0x0f) == ECT_MBXT_SOE) &&
                   (aSoEp->opCode == ECT_SOE_WRITERES) &&
                   (aSoEp->error == 0) &&
                   (aSoEp->driveNo == driveNo) &&
                   (aSoEp->elementflags == elementflags))
               {
                  /* SoE write succeeded */
               }
               /* other slave response */
               else
               {
                  if (((aSoEp->MbxHeader.mbxtype & 0x0f) == ECT_MBXT_SOE) &&
                      (aSoEp->opCode == ECT_SOE_READRES) &&
                      (aSoEp->error == 1))
                  {
                     mp = (uint8 *)&MbxIn + (etohs(aSoEp->MbxHeader.length) + sizeof(ec_mbxheadert) - sizeof(uint16));
                     errorcode = (uint16 *)mp;
                     ecx_SoEerror(context, slave, idn, *errorcode);
                  }
                  else
                  {
                     ecx_packeterror(context, slave, idn, 0, 1); /* Unexpected frame returned */
                  }
                  wkc = 0;
               }
            }
            else
            {
               ecx_packeterror(context, slave, idn, 0, 4); /* no response */
            }
         }
      }
   }
   return wkc;
}

/** SoE read AT and MTD mapping.
 *
 * SoE has standard indexes defined for mapping. This function
 * tries to read them and collect a full input and output mapping size
 * of designated slave.
 *
 * @param[in]  context = context struct
 * @param[in]  slave   = Slave number
 * @param[out] Osize   = Size in bits of output mapping (MTD) found
 * @param[out] Isize   = Size in bits of input mapping (AT) found
 * @return >0 if mapping succesful.
 */
int ecx_readIDNmap(ecx_contextt *context, uint16 slave, int *Osize, int *Isize)
{
   int retVal = 0;
   int   wkc;
   int psize;
   int driveNr;
   uint16 entries, itemcount;
   ec_SoEmappingt     SoEmapping;
   ec_SoEattributet   SoEattribute;

   *Isize = 0;
   *Osize = 0;
   for(driveNr = 0; driveNr < EC_SOE_MAX_DRIVES; driveNr++)
   {
      psize = sizeof(SoEmapping);
      /* read output mapping via SoE */
      wkc = ecx_SoEread(context, slave, driveNr, EC_SOE_VALUE_B, EC_IDN_MDTCONFIG, &psize, &SoEmapping, EC_TIMEOUTRXM);
      if ((wkc > 0) && (psize >= 4) && ((entries = etohs(SoEmapping.currentlength) / 2) > 0) && (entries <= EC_SOE_MAXMAPPING))
      {
         /* command word (uint16) is always mapped but not in list */
         *Osize = 16;
         for (itemcount = 0 ; itemcount < entries ; itemcount++)
         {
            psize = sizeof(SoEattribute);
            /* read attribute of each IDN in mapping list */
            wkc = ecx_SoEread(context, slave, driveNr, EC_SOE_ATTRIBUTE_B, SoEmapping.idn[itemcount], &psize, &SoEattribute, EC_TIMEOUTRXM);
            if ((wkc > 0) && (!SoEattribute.list))
            {
               /* length : 0 = 8bit, 1 = 16bit .... */
               *Osize += (int)8 << SoEattribute.length;
            }
         }
      }
      psize = sizeof(SoEmapping);
      /* read input mapping via SoE */
      wkc = ecx_SoEread(context, slave, driveNr, EC_SOE_VALUE_B, EC_IDN_ATCONFIG, &psize, &SoEmapping, EC_TIMEOUTRXM);
      if ((wkc > 0) && (psize >= 4) && ((entries = etohs(SoEmapping.currentlength) / 2) > 0) && (entries <= EC_SOE_MAXMAPPING))
      {
         /* status word (uint16) is always mapped but not in list */
         *Isize = 16;
         for (itemcount = 0 ; itemcount < entries ; itemcount++)
         {
            psize = sizeof(SoEattribute);
            /* read attribute of each IDN in mapping list */
            wkc = ecx_SoEread(context, slave, driveNr, EC_SOE_ATTRIBUTE_B, SoEmapping.idn[itemcount], &psize, &SoEattribute, EC_TIMEOUTRXM);
            if ((wkc > 0) && (!SoEattribute.list))
            {
               /* length : 0 = 8bit, 1 = 16bit .... */
               *Isize += (int)8 << SoEattribute.length;
            }
         }
      }
   }

   /* found some I/O bits ? */
   if ((*Isize > 0) || (*Osize > 0))
   {
      retVal = 1;
   }
   return retVal;
}

#ifdef EC_VER1
int ec_SoEread(uint16 slave, uint8 driveNo, uint8 elementflags, uint16 idn, int *psize, void *p, int timeout)
{
   return ecx_SoEread(&ecx_context, slave, driveNo, elementflags, idn, psize, p, timeout);
}

int ec_SoEwrite(uint16 slave, uint8 driveNo, uint8 elementflags, uint16 idn, int psize, void *p, int timeout)
{
   return ecx_SoEwrite(&ecx_context, slave, driveNo, elementflags, idn, psize, p, timeout);
}

int ec_readIDNmap(uint16 slave, int *Osize, int *Isize)
{
   return ecx_readIDNmap(&ecx_context, slave, Osize, Isize);
}
#endif
