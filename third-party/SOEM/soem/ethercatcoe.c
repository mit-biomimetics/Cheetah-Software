/*
 * Simple Open EtherCAT Master Library
 *
 * File    : ethercatcoe.c
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
 * CAN over EtherCAT (CoE) module.
 *
 * SDO read / write and SDO service functions
 */

#include <stdio.h>
#include <string.h>
#include "osal.h"
#include "oshw.h"
#include "ethercattype.h"
#include "ethercatbase.h"
#include "ethercatmain.h"
#include "ethercatcoe.h"

/** SDO structure, not to be confused with EcSDOserviceT */
PACKED_BEGIN
typedef struct PACKED
{
   ec_mbxheadert   MbxHeader;
   uint16          CANOpen;
   uint8           Command;
   uint16          Index;
   uint8           SubIndex;
   union
   {
      uint8   bdata[0x200]; /* variants for easy data access */
      uint16  wdata[0x100];
      uint32  ldata[0x80];
   };
} ec_SDOt;
PACKED_END

/** SDO service structure */
PACKED_BEGIN
typedef struct PACKED
{
   ec_mbxheadert   MbxHeader;
   uint16          CANOpen;
   uint8           Opcode;
   uint8           Reserved;
   uint16          Fragments;
   union
   {
      uint8   bdata[0x200]; /* variants for easy data access */
      uint16  wdata[0x100];
      uint32  ldata[0x80];
   };
} ec_SDOservicet;
PACKED_END

/** Report SDO error.
 *
 * @param[in]  context    = context struct
 * @param[in]  Slave      = Slave number
 * @param[in]  Index      = Index that generated error
 * @param[in]  SubIdx     = Subindex that generated error
 * @param[in]  AbortCode  = Abortcode, see EtherCAT documentation for list
 */
void ecx_SDOerror(ecx_contextt *context, uint16 Slave, uint16 Index, uint8 SubIdx, int32 AbortCode)
{
   ec_errort Ec;

   memset(&Ec, 0, sizeof(Ec));
   Ec.Time = osal_current_time();
   Ec.Slave = Slave;
   Ec.Index = Index;
   Ec.SubIdx = SubIdx;
   *(context->ecaterror) = TRUE;
   Ec.Etype = EC_ERR_TYPE_SDO_ERROR;
   Ec.AbortCode = AbortCode;
   ecx_pusherror(context, &Ec);
}

/** Report SDO info error
 *
 * @param[in]  context    = context struct
 * @param[in]  Slave      = Slave number
 * @param[in]  Index      = Index that generated error
 * @param[in]  SubIdx     = Subindex that generated error
 * @param[in]  AbortCode  = Abortcode, see EtherCAT documentation for list
 */
static void ecx_SDOinfoerror(ecx_contextt *context, uint16 Slave, uint16 Index, uint8 SubIdx, int32 AbortCode)
{
   ec_errort Ec;

   memset(&Ec, 0, sizeof(Ec));
   Ec.Slave = Slave;
   Ec.Index = Index;
   Ec.SubIdx = SubIdx;
   *(context->ecaterror) = TRUE;
   Ec.Etype = EC_ERR_TYPE_SDOINFO_ERROR;
   Ec.AbortCode = AbortCode;
   ecx_pusherror(context, &Ec);
}

/** CoE SDO read, blocking. Single subindex or Complete Access.
 *
 * Only a "normal" upload request is issued. If the requested parameter is <= 4bytes
 * then a "expedited" response is returned, otherwise a "normal" response. If a "normal"
 * response is larger than the mailbox size then the response is segmented. The function
 * will combine all segments and copy them to the parameter buffer.
 *
 * @param[in]  context    = context struct
 * @param[in]  slave      = Slave number
 * @param[in]  index      = Index to read
 * @param[in]  subindex   = Subindex to read, must be 0 or 1 if CA is used.
 * @param[in]  CA         = FALSE = single subindex. TRUE = Complete Access, all subindexes read.
 * @param[in,out] psize   = Size in bytes of parameter buffer, returns bytes read from SDO.
 * @param[out] p          = Pointer to parameter buffer
 * @param[in]  timeout    = Timeout in us, standard is EC_TIMEOUTRXM
 * @return Workcounter from last slave response
 */
int ecx_SDOread(ecx_contextt *context, uint16 slave, uint16 index, uint8 subindex,
               boolean CA, int *psize, void *p, int timeout)
{
   ec_SDOt *SDOp, *aSDOp;
   uint16 bytesize, Framedatasize;
   int wkc;
   int32 SDOlen;
   uint8 *bp;
   uint8 *hp;
   ec_mbxbuft MbxIn, MbxOut;
   uint8 cnt, toggle;
   boolean NotLast;

   ec_clearmbx(&MbxIn);
   /* Empty slave out mailbox if something is in. Timout set to 0 */
   wkc = ecx_mbxreceive(context, slave, (ec_mbxbuft *)&MbxIn, 0);
   ec_clearmbx(&MbxOut);
   aSDOp = (ec_SDOt *)&MbxIn;
   SDOp = (ec_SDOt *)&MbxOut;
   SDOp->MbxHeader.length = htoes(0x000a);
   SDOp->MbxHeader.address = htoes(0x0000);
   SDOp->MbxHeader.priority = 0x00;
   /* get new mailbox count value, used as session handle */
   cnt = ec_nextmbxcnt(context->slavelist[slave].mbx_cnt);
   context->slavelist[slave].mbx_cnt = cnt;
   SDOp->MbxHeader.mbxtype = ECT_MBXT_COE + (cnt << 4); /* CoE */
   SDOp->CANOpen = htoes(0x000 + (ECT_COES_SDOREQ << 12)); /* number 9bits service upper 4 bits (SDO request) */
   if (CA)
   {
      SDOp->Command = ECT_SDO_UP_REQ_CA; /* upload request complete access */
   }
   else
   {
      SDOp->Command = ECT_SDO_UP_REQ; /* upload request normal */
   }
   SDOp->Index = htoes(index);
   if (CA && (subindex > 1))
   {
      subindex = 1;
   }
   SDOp->SubIndex = subindex;
   SDOp->ldata[0] = 0;
   /* send CoE request to slave */
   wkc = ecx_mbxsend(context, slave, (ec_mbxbuft *)&MbxOut, EC_TIMEOUTTXM);
   if (wkc > 0) /* succeeded to place mailbox in slave ? */
   {
      /* clean mailboxbuffer */
      ec_clearmbx(&MbxIn);
      /* read slave response */
      wkc = ecx_mbxreceive(context, slave, (ec_mbxbuft *)&MbxIn, timeout);
      if (wkc > 0) /* succeeded to read slave response ? */
      {
         /* slave response should be CoE, SDO response and the correct index */
         if (((aSDOp->MbxHeader.mbxtype & 0x0f) == ECT_MBXT_COE) &&
             ((etohs(aSDOp->CANOpen) >> 12) == ECT_COES_SDORES) &&
              (aSDOp->Index == SDOp->Index))
         {
            if ((aSDOp->Command & 0x02) > 0)
            {
               /* expedited frame response */
               bytesize = 4 - ((aSDOp->Command >> 2) & 0x03);
               if (*psize >= bytesize) /* parameter buffer big enough ? */
               {
                  /* copy parameter in parameter buffer */
                  memcpy(p, &aSDOp->ldata[0], bytesize);
                  /* return the real parameter size */
                  *psize = bytesize;
               }
               else
               {
                  wkc = 0;
                  ecx_packeterror(context, slave, index, subindex, 3); /*  data container too small for type */
               }
            }
            else
            { /* normal frame response */
               SDOlen = etohl(aSDOp->ldata[0]);
               /* Does parameter fit in parameter buffer ? */
               if (SDOlen <= *psize)
               {
                  bp = p;
                  hp = p;
                  /* calculate mailbox transfer size */
                  Framedatasize = (etohs(aSDOp->MbxHeader.length) - 10);
                  if (Framedatasize < SDOlen) /* transfer in segments? */
                  {
                     /* copy parameter data in parameter buffer */
                     memcpy(hp, &aSDOp->ldata[1], Framedatasize);
                     /* increment buffer pointer */
                     hp += Framedatasize;
                     *psize = Framedatasize;
                     NotLast = TRUE;
                     toggle= 0x00;
                     while (NotLast) /* segmented transfer */
                     {
                        SDOp = (ec_SDOt *)&MbxOut;
                        SDOp->MbxHeader.length = htoes(0x000a);
                        SDOp->MbxHeader.address = htoes(0x0000);
                        SDOp->MbxHeader.priority = 0x00;
                        cnt = ec_nextmbxcnt(context->slavelist[slave].mbx_cnt);
                        context->slavelist[slave].mbx_cnt = cnt;
                        SDOp->MbxHeader.mbxtype = ECT_MBXT_COE + (cnt << 4); /* CoE */
                        SDOp->CANOpen = htoes(0x000 + (ECT_COES_SDOREQ << 12)); /* number 9bits service upper 4 bits (SDO request) */
                        SDOp->Command = ECT_SDO_SEG_UP_REQ + toggle; /* segment upload request */
                        SDOp->Index = htoes(index);
                        SDOp->SubIndex = subindex;
                        SDOp->ldata[0] = 0;
                        /* send segmented upload request to slave */
                        wkc = ecx_mbxsend(context, slave, (ec_mbxbuft *)&MbxOut, EC_TIMEOUTTXM);
                        /* is mailbox transfered to slave ? */
                        if (wkc > 0)
                        {
                           ec_clearmbx(&MbxIn);
                           /* read slave response */
                           wkc = ecx_mbxreceive(context, slave, (ec_mbxbuft *)&MbxIn, timeout);
                           /* has slave responded ? */
                           if (wkc > 0)
                           {
                              /* slave response should be CoE, SDO response */
                              if ((((aSDOp->MbxHeader.mbxtype & 0x0f) == ECT_MBXT_COE) &&
                                   ((etohs(aSDOp->CANOpen) >> 12) == ECT_COES_SDORES) &&
                                   ((aSDOp->Command & 0xe0) == 0x00)))
                                        {
                                 /* calculate mailbox transfer size */
                                 Framedatasize = etohs(aSDOp->MbxHeader.length) - 3;
                                 if ((aSDOp->Command & 0x01) > 0)
                                 { /* last segment */
                                    NotLast = FALSE;
                                    if (Framedatasize == 7)
                                       /* substract unused bytes from frame */
                                       Framedatasize = Framedatasize - ((aSDOp->Command & 0x0e) >> 1);
                                    /* copy to parameter buffer */
                                    memcpy(hp, &(aSDOp->Index), Framedatasize);
                                 }
                                 else /* segments follow */
                                 {
                                    /* copy to parameter buffer */
                                    memcpy(hp, &(aSDOp->Index), Framedatasize);
                                    /* increment buffer pointer */
                                    hp += Framedatasize;
                                 }
                                 /* update parametersize */
                                 *psize += Framedatasize;
                              }
                              /* unexpected frame returned from slave */
                              else
                              {
                                 NotLast = FALSE;
                                 if ((aSDOp->Command) == ECT_SDO_ABORT) /* SDO abort frame received */
                                    ecx_SDOerror(context, slave, index, subindex, etohl(aSDOp->ldata[0]));
                                 else
                                    ecx_packeterror(context, slave, index, subindex, 1); /* Unexpected frame returned */
                                 wkc = 0;
                              }
                           }
                        }
                        toggle = toggle ^ 0x10; /* toggle bit for segment request */
                     }
                  }
                  /* non segmented transfer */
                  else
                  {
                     /* copy to parameter buffer */
                     memcpy(bp, &aSDOp->ldata[1], SDOlen);
                     *psize = SDOlen;
                  }
               }
               /* parameter buffer too small */
               else
               {
                  wkc = 0;
                  ecx_packeterror(context, slave, index, subindex, 3); /*  data container too small for type */
               }
            }
         }
         /* other slave response */
         else
         {
            if ((aSDOp->Command) == ECT_SDO_ABORT) /* SDO abort frame received */
            {
               ecx_SDOerror(context, slave, index, subindex, etohl(aSDOp->ldata[0]));
            }
            else
            {
               ecx_packeterror(context, slave, index, subindex, 1); /* Unexpected frame returned */
            }
            wkc = 0;
         }
      }
   }
   return wkc;
}

/** CoE SDO write, blocking. Single subindex or Complete Access.
 *
 * A "normal" download request is issued, unless we have
 * small data, then a "expedited" transfer is used. If the parameter is larger than
 * the mailbox size then the download is segmented. The function will split the
 * parameter data in segments and send them to the slave one by one.
 *
 * @param[in]  context    = context struct
 * @param[in]  Slave      = Slave number
 * @param[in]  Index      = Index to write
 * @param[in]  SubIndex   = Subindex to write, must be 0 or 1 if CA is used.
 * @param[in]  CA         = FALSE = single subindex. TRUE = Complete Access, all subindexes written.
 * @param[in]  psize      = Size in bytes of parameter buffer.
 * @param[out] p          = Pointer to parameter buffer
 * @param[in]  Timeout    = Timeout in us, standard is EC_TIMEOUTRXM
 * @return Workcounter from last slave response
 */
int ecx_SDOwrite(ecx_contextt *context, uint16 Slave, uint16 Index, uint8 SubIndex,
                boolean CA, int psize, void *p, int Timeout)
{
   ec_SDOt *SDOp, *aSDOp;
   int wkc, maxdata;
   ec_mbxbuft MbxIn, MbxOut;
   uint8 cnt, toggle;
   uint16 framedatasize;
   boolean  NotLast;
   uint8 *hp;

   ec_clearmbx(&MbxIn);
   /* Empty slave out mailbox if something is in. Timout set to 0 */
   wkc = ecx_mbxreceive(context, Slave, (ec_mbxbuft *)&MbxIn, 0);
   ec_clearmbx(&MbxOut);
   aSDOp = (ec_SDOt *)&MbxIn;
   SDOp = (ec_SDOt *)&MbxOut;
   maxdata = context->slavelist[Slave].mbx_l - 0x10; /* data section=mailbox size - 6 mbx - 2 CoE - 8 sdo req */
   /* if small data use expedited transfer */
   if ((psize <= 4) && !CA)
   {
      SDOp->MbxHeader.length = htoes(0x000a);
      SDOp->MbxHeader.address = htoes(0x0000);
      SDOp->MbxHeader.priority = 0x00;
      /* get new mailbox counter, used for session handle */
      cnt = ec_nextmbxcnt(context->slavelist[Slave].mbx_cnt);
      context->slavelist[Slave].mbx_cnt = cnt;
      SDOp->MbxHeader.mbxtype = ECT_MBXT_COE + (cnt << 4); /* CoE */
      SDOp->CANOpen = htoes(0x000 + (ECT_COES_SDOREQ << 12)); /* number 9bits service upper 4 bits */
      SDOp->Command = ECT_SDO_DOWN_EXP | (((4 - psize) << 2) & 0x0c); /* expedited SDO download transfer */
      SDOp->Index = htoes(Index);
      SDOp->SubIndex = SubIndex;
      hp = p;
      /* copy parameter data to mailbox */
      memcpy(&SDOp->ldata[0], hp, psize);
      /* send mailbox SDO download request to slave */
      wkc = ecx_mbxsend(context, Slave, (ec_mbxbuft *)&MbxOut, EC_TIMEOUTTXM);
      if (wkc > 0)
      {
         ec_clearmbx(&MbxIn);
         /* read slave response */
         wkc = ecx_mbxreceive(context, Slave, (ec_mbxbuft *)&MbxIn, Timeout);
         if (wkc > 0)
         {
            /* response should be CoE, SDO response, correct index and subindex */
            if (((aSDOp->MbxHeader.mbxtype & 0x0f) == ECT_MBXT_COE) &&
                ((etohs(aSDOp->CANOpen) >> 12) == ECT_COES_SDORES) &&
                 (aSDOp->Index == SDOp->Index) &&
                 (aSDOp->SubIndex == SDOp->SubIndex))
            {
                 /* all OK */
            }
            /* unexpected response from slave */
            else
            {
               if (aSDOp->Command == ECT_SDO_ABORT) /* SDO abort frame received */
               {
                  ecx_SDOerror(context, Slave, Index, SubIndex, etohl(aSDOp->ldata[0]));
               }
               else
               {
                  ecx_packeterror(context, Slave, Index, SubIndex, 1); /* Unexpected frame returned */
               }
               wkc = 0;
            }
         }
      }
   }
   else
   {
      framedatasize = psize;
      NotLast = FALSE;
      if (framedatasize > maxdata)
      {
         framedatasize = maxdata;  /*  segmented transfer needed  */
         NotLast = TRUE;
      }
      SDOp->MbxHeader.length = htoes(0x0a + framedatasize);
      SDOp->MbxHeader.address = htoes(0x0000);
      SDOp->MbxHeader.priority = 0x00;
      /* get new mailbox counter, used for session handle */
      cnt = ec_nextmbxcnt(context->slavelist[Slave].mbx_cnt);
      context->slavelist[Slave].mbx_cnt = cnt;
      SDOp->MbxHeader.mbxtype = ECT_MBXT_COE + (cnt << 4); /* CoE */
      SDOp->CANOpen = htoes(0x000 + (ECT_COES_SDOREQ << 12)); /* number 9bits service upper 4 bits */
      if (CA)
      {
         SDOp->Command = ECT_SDO_DOWN_INIT_CA; /* Complete Access, normal SDO init download transfer */
      }
      else
      {
         SDOp->Command = ECT_SDO_DOWN_INIT; /* normal SDO init download transfer */
      }
      SDOp->Index = htoes(Index);
      SDOp->SubIndex = SubIndex;
      if (CA && (SubIndex > 1))
      {
         SDOp->SubIndex = 1;
      }
      SDOp->ldata[0] = htoel(psize);
      hp = p;
      /* copy parameter data to mailbox */
      memcpy(&SDOp->ldata[1], hp, framedatasize);
      hp += framedatasize;
      psize -= framedatasize;
      /* send mailbox SDO download request to slave */
      wkc = ecx_mbxsend(context, Slave, (ec_mbxbuft *)&MbxOut, EC_TIMEOUTTXM);
      if (wkc > 0)
      {
         ec_clearmbx(&MbxIn);
         /* read slave response */
         wkc = ecx_mbxreceive(context, Slave, (ec_mbxbuft *)&MbxIn, Timeout);
         if (wkc > 0)
         {
            /* response should be CoE, SDO response, correct index and subindex */
            if (((aSDOp->MbxHeader.mbxtype & 0x0f) == ECT_MBXT_COE) &&
                ((etohs(aSDOp->CANOpen) >> 12) == ECT_COES_SDORES) &&
                 (aSDOp->Index == SDOp->Index) &&
                 (aSDOp->SubIndex == SDOp->SubIndex))
            {
               /* all ok */
               maxdata += 7;
               toggle = 0;
               /* repeat while segments left */
               while (NotLast)
               {
                  SDOp = (ec_SDOt *)&MbxOut;
                  framedatasize = psize;
                  NotLast = FALSE;
                  SDOp->Command = 0x01; /* last segment */
                  if (framedatasize > maxdata)
                  {
                     framedatasize = maxdata;  /*  more segments needed  */
                     NotLast = TRUE;
                     SDOp->Command = 0x00; /* segments follow */
                  }
                  if (!NotLast && (framedatasize < 7))
                  {
                     SDOp->MbxHeader.length = htoes(0x0a); /* minimum size */
                     SDOp->Command = 0x01 + ((7 - framedatasize) << 1); /* last segment reduced octets */
                  }
                  else
                  {
                     SDOp->MbxHeader.length = htoes(framedatasize + 3); /* data + 2 CoE + 1 SDO */
                  }
                  SDOp->MbxHeader.address = htoes(0x0000);
                  SDOp->MbxHeader.priority = 0x00;
                  /* get new mailbox counter value */
                  cnt = ec_nextmbxcnt(context->slavelist[Slave].mbx_cnt);
                  context->slavelist[Slave].mbx_cnt = cnt;
                  SDOp->MbxHeader.mbxtype = ECT_MBXT_COE + (cnt << 4); /* CoE */
                  SDOp->CANOpen = htoes(0x000 + (ECT_COES_SDOREQ << 12)); /* number 9bits service upper 4 bits (SDO request) */
                  SDOp->Command = SDOp->Command + toggle; /* add toggle bit to command byte */
                  /* copy parameter data to mailbox */
                  memcpy(&SDOp->Index, hp, framedatasize);
                  /* update parameter buffer pointer */
                  hp += framedatasize;
                  psize -= framedatasize;
                  /* send SDO download request */
                  wkc = ecx_mbxsend(context, Slave, (ec_mbxbuft *)&MbxOut, EC_TIMEOUTTXM);
                  if (wkc > 0)
                  {
                     ec_clearmbx(&MbxIn);
                     /* read slave response */
                     wkc = ecx_mbxreceive(context, Slave, (ec_mbxbuft *)&MbxIn, Timeout);
                     if (wkc > 0)
                     {
                        if (((aSDOp->MbxHeader.mbxtype & 0x0f) == ECT_MBXT_COE) &&
                            ((etohs(aSDOp->CANOpen) >> 12) == ECT_COES_SDORES) &&
                            ((aSDOp->Command & 0xe0) == 0x20))
                        {
                                   /* all OK, nothing to do */
                        }
                        else
                        {
                           if (aSDOp->Command == ECT_SDO_ABORT) /* SDO abort frame received */
                           {
                              ecx_SDOerror(context, Slave, Index, SubIndex, etohl(aSDOp->ldata[0]));
                           }
                           else
                           {
                              ecx_packeterror(context, Slave, Index, SubIndex, 1); /* Unexpected frame returned */
                           }
                           wkc = 0;
                           NotLast = FALSE;
                        }
                     }
                  }
                  toggle = toggle ^ 0x10; /* toggle bit for segment request */
               }
            }
            /* unexpected response from slave */
            else
            {
               if (aSDOp->Command == ECT_SDO_ABORT) /* SDO abort frame received */
               {
                  ecx_SDOerror(context, Slave, Index, SubIndex, etohl(aSDOp->ldata[0]));
               }
               else
               {
                  ecx_packeterror(context, Slave, Index, SubIndex, 1); /* Unexpected frame returned */
               }
               wkc = 0;
            }
         }
      }
   }

   return wkc;
}

/** CoE RxPDO write, blocking.
 *
 * A RxPDO download request is issued.
 *
 * @param[in]  context       = context struct
 * @param[in]  Slave         = Slave number
 * @param[in]  RxPDOnumber   = Related RxPDO number
 * @param[in]  psize         = Size in bytes of PDO buffer.
 * @param[out] p             = Pointer to PDO buffer
 * @return Workcounter from last slave response
 */
int ecx_RxPDO(ecx_contextt *context, uint16 Slave, uint16 RxPDOnumber, int psize, void *p)
{
   ec_SDOt *SDOp;
   int wkc, maxdata;
   ec_mbxbuft MbxIn, MbxOut;
   uint8 cnt;
   uint16 framedatasize;

   ec_clearmbx(&MbxIn);
   /* Empty slave out mailbox if something is in. Timout set to 0 */
   wkc = ecx_mbxreceive(context, Slave, (ec_mbxbuft *)&MbxIn, 0);
   ec_clearmbx(&MbxOut);
   SDOp = (ec_SDOt *)&MbxOut;
   maxdata = context->slavelist[Slave].mbx_l - 0x08; /* data section=mailbox size - 6 mbx - 2 CoE */
   framedatasize = psize;
   if (framedatasize > maxdata)
   {
      framedatasize = maxdata;  /*  limit transfer */
   }
   SDOp->MbxHeader.length = htoes(0x02 + framedatasize);
   SDOp->MbxHeader.address = htoes(0x0000);
   SDOp->MbxHeader.priority = 0x00;
   /* get new mailbox counter, used for session handle */
   cnt = ec_nextmbxcnt(context->slavelist[Slave].mbx_cnt);
   context->slavelist[Slave].mbx_cnt = cnt;
   SDOp->MbxHeader.mbxtype = ECT_MBXT_COE + (cnt << 4); /* CoE */
   SDOp->CANOpen = htoes((RxPDOnumber & 0x01ff) + (ECT_COES_RXPDO << 12)); /* number 9bits service upper 4 bits */
   /* copy PDO data to mailbox */
   memcpy(&SDOp->Command, p, framedatasize);
   /* send mailbox RxPDO request to slave */
   wkc = ecx_mbxsend(context, Slave, (ec_mbxbuft *)&MbxOut, EC_TIMEOUTTXM);

   return wkc;
}

/** CoE TxPDO read remote request, blocking.
 *
 * A RxPDO download request is issued.
 *
 * @param[in]  context       = context struct
 * @param[in]  slave         = Slave number
 * @param[in]  TxPDOnumber   = Related TxPDO number
 * @param[in,out] psize      = Size in bytes of PDO buffer, returns bytes read from PDO.
 * @param[out] p             = Pointer to PDO buffer
 * @param[in]  timeout       = Timeout in us, standard is EC_TIMEOUTRXM
 * @return Workcounter from last slave response
 */
int ecx_TxPDO(ecx_contextt *context, uint16 slave, uint16 TxPDOnumber , int *psize, void *p, int timeout)
{
   ec_SDOt *SDOp, *aSDOp;
   int wkc;
   ec_mbxbuft MbxIn, MbxOut;
   uint8 cnt;
   uint16 framedatasize;

   ec_clearmbx(&MbxIn);
   /* Empty slave out mailbox if something is in. Timout set to 0 */
   wkc = ecx_mbxreceive(context, slave, (ec_mbxbuft *)&MbxIn, 0);
   ec_clearmbx(&MbxOut);
   aSDOp = (ec_SDOt *)&MbxIn;
   SDOp = (ec_SDOt *)&MbxOut;
   SDOp->MbxHeader.length = htoes(0x02);
   SDOp->MbxHeader.address = htoes(0x0000);
   SDOp->MbxHeader.priority = 0x00;
   /* get new mailbox counter, used for session handle */
   cnt = ec_nextmbxcnt(context->slavelist[slave].mbx_cnt);
   context->slavelist[slave].mbx_cnt = cnt;
   SDOp->MbxHeader.mbxtype = ECT_MBXT_COE + (cnt << 4); /* CoE */
   SDOp->CANOpen = htoes((TxPDOnumber & 0x01ff) + (ECT_COES_TXPDO_RR << 12)); /* number 9bits service upper 4 bits */
   wkc = ecx_mbxsend(context, slave, (ec_mbxbuft *)&MbxOut, EC_TIMEOUTTXM);
   if (wkc > 0)
   {
      /* clean mailboxbuffer */
      ec_clearmbx(&MbxIn);
      /* read slave response */
      wkc = ecx_mbxreceive(context, slave, (ec_mbxbuft *)&MbxIn, timeout);
      if (wkc > 0) /* succeeded to read slave response ? */
      {
         /* slave response should be CoE, TxPDO */
         if (((aSDOp->MbxHeader.mbxtype & 0x0f) == ECT_MBXT_COE) &&
             ((etohs(aSDOp->CANOpen) >> 12) == ECT_COES_TXPDO))
         {
            /* TxPDO response */
            framedatasize = (aSDOp->MbxHeader.length - 2);
            if (*psize >= framedatasize) /* parameter buffer big enough ? */
            {
               /* copy parameter in parameter buffer */
               memcpy(p, &aSDOp->Command, framedatasize);
               /* return the real parameter size */
               *psize = framedatasize;
            }
            /* parameter buffer too small */
            else
            {
               wkc = 0;
               ecx_packeterror(context, slave, 0, 0, 3); /*  data container too small for type */
            }
         }
         /* other slave response */
         else
         {
            if ((aSDOp->Command) == ECT_SDO_ABORT) /* SDO abort frame received */
            {
               ecx_SDOerror(context, slave, 0, 0, etohl(aSDOp->ldata[0]));
            }
            else
            {
               ecx_packeterror(context, slave, 0, 0, 1); /* Unexpected frame returned */
            }
            wkc = 0;
         }
      }
   }

   return wkc;
}

/** Read PDO assign structure
 * @param[in]  context       = context struct
 * @param[in]  Slave         = Slave number
 * @param[in]  PDOassign     = PDO assign object
 * @return total bitlength of PDO assign
 */
int ecx_readPDOassign(ecx_contextt *context, uint16 Slave, uint16 PDOassign)
{
   uint16 idxloop, nidx, subidxloop, rdat, idx, subidx;
   uint8 subcnt;
   int wkc, bsize = 0, rdl;
   int32 rdat2;

   rdl = sizeof(rdat); rdat = 0;
   /* read PDO assign subindex 0 ( = number of PDO's) */
   wkc = ecx_SDOread(context, Slave, PDOassign, 0x00, FALSE, &rdl, &rdat, EC_TIMEOUTRXM);
   rdat = etohs(rdat);
   /* positive result from slave ? */
   if ((wkc > 0) && (rdat > 0))
   {
      /* number of available sub indexes */
      nidx = rdat;
      bsize = 0;
      /* read all PDO's */
      for (idxloop = 1; idxloop <= nidx; idxloop++)
      {
         rdl = sizeof(rdat); rdat = 0;
         /* read PDO assign */
         wkc = ecx_SDOread(context, Slave, PDOassign, (uint8)idxloop, FALSE, &rdl, &rdat, EC_TIMEOUTRXM);
         /* result is index of PDO */
         idx = etohl(rdat);
         if (idx > 0)
         {
            rdl = sizeof(subcnt); subcnt = 0;
            /* read number of subindexes of PDO */
            wkc = ecx_SDOread(context, Slave,idx, 0x00, FALSE, &rdl, &subcnt, EC_TIMEOUTRXM);
            subidx = subcnt;
            /* for each subindex */
            for (subidxloop = 1; subidxloop <= subidx; subidxloop++)
            {
               rdl = sizeof(rdat2); rdat2 = 0;
               /* read SDO that is mapped in PDO */
               wkc = ecx_SDOread(context, Slave, idx, (uint8)subidxloop, FALSE, &rdl, &rdat2, EC_TIMEOUTRXM);
               rdat2 = etohl(rdat2);
               /* extract bitlength of SDO */
               if (LO_BYTE(rdat2) < 0xff)
               {
                  bsize += LO_BYTE(rdat2);
               }
               else
               {
                  rdl = sizeof(rdat); rdat = htoes(0xff);
                  /* read Object Entry in Object database */
//                  wkc = ec_readOEsingle(idx, (uint8)SubCount, pODlist, pOElist);
                  bsize += etohs(rdat);
               }
            }
         }
      }
   }
   /* return total found bitlength (PDO) */
   return bsize;
}

/** Read PDO assign structure in Complete Access mode
 * @param[in]  context       = context struct
 * @param[in]  Slave         = Slave number
 * @param[in]  Thread_n      = Calling thread index
 * @param[in]  PDOassign     = PDO assign object
 * @return total bitlength of PDO assign
 */
int ecx_readPDOassignCA(ecx_contextt *context, uint16 Slave, int Thread_n,
      uint16 PDOassign)
{
   uint16 idxloop, nidx, subidxloop, idx, subidx;
   int wkc, bsize = 0, rdl;

   /* find maximum size of PDOassign buffer */
   rdl = sizeof(ec_PDOassignt);
   context->PDOassign[Thread_n].n=0;
   /* read rxPDOassign in CA mode, all subindexes are read in one struct */
   wkc = ecx_SDOread(context, Slave, PDOassign, 0x00, TRUE, &rdl,
         &(context->PDOassign[Thread_n]), EC_TIMEOUTRXM);
   /* positive result from slave ? */
   if ((wkc > 0) && (context->PDOassign[Thread_n].n > 0))
   {
      nidx = context->PDOassign[Thread_n].n;
      bsize = 0;
      /* for each PDO do */
      for (idxloop = 1; idxloop <= nidx; idxloop++)
      {
         /* get index from PDOassign struct */
         idx = etohs(context->PDOassign[Thread_n].index[idxloop - 1]);
         if (idx > 0)
         {
            rdl = sizeof(ec_PDOdesct); context->PDOdesc[Thread_n].n = 0;
            /* read SDO's that are mapped in PDO, CA mode */
            wkc = ecx_SDOread(context, Slave,idx, 0x00, TRUE, &rdl,
                  &(context->PDOdesc[Thread_n]), EC_TIMEOUTRXM);
            subidx = context->PDOdesc[Thread_n].n;
            /* extract all bitlengths of SDO's */
            for (subidxloop = 1; subidxloop <= subidx; subidxloop++)
            {
               bsize += LO_BYTE(etohl(context->PDOdesc[Thread_n].PDO[subidxloop -1]));
            }
         }
      }
   }

   /* return total found bitlength (PDO) */
   return bsize;
}

/** CoE read PDO mapping.
 *
 * CANopen has standard indexes defined for PDO mapping. This function
 * tries to read them and collect a full input and output mapping size
 * of designated slave.
 *
 * Principal structure in slave:\n
 * 1C00:00 is number of SM defined\n
 * 1C00:01 SM0 type -> 1C10\n
 * 1C00:02 SM1 type -> 1C11\n
 * 1C00:03 SM2 type -> 1C12\n
 * 1C00:04 SM3 type -> 1C13\n
 * Type 0 = unused, 1 = mailbox in, 2 = mailbox out,
 * 3 = outputs (RxPDO), 4 = inputs (TxPDO).
 *
 * 1C12:00 is number of PDO's defined for SM2\n
 * 1C12:01 PDO assign SDO #1 -> f.e. 1A00\n
 * 1C12:02 PDO assign SDO #2 -> f.e. 1A04\
 *
 * 1A00:00 is number of object defined for this PDO\n
 * 1A00:01 object mapping #1, f.e. 60100710 (SDO 6010 SI 07 bitlength 0x10)
 *
 * @param[in]  context = context struct
 * @param[in]  Slave   = Slave number
 * @param[out] Osize   = Size in bits of output mapping (rxPDO) found
 * @param[out] Isize   = Size in bits of input mapping (txPDO) found
 * @return >0 if mapping succesful.
 */
int ecx_readPDOmap(ecx_contextt *context, uint16 Slave, int *Osize, int *Isize)
{
   int wkc, rdl;
   int retVal = 0;
   uint8 nSM, iSM, tSM;
   int Tsize;
   uint8 SMt_bug_add;

   *Isize = 0;
   *Osize = 0;
   SMt_bug_add = 0;
   rdl = sizeof(nSM); nSM = 0;
   /* read SyncManager Communication Type object count */
   wkc = ecx_SDOread(context, Slave, ECT_SDO_SMCOMMTYPE, 0x00, FALSE, &rdl, &nSM, EC_TIMEOUTRXM);
   /* positive result from slave ? */
   if ((wkc > 0) && (nSM > 2))
   {
      /* limit to maximum number of SM defined, if true the slave can't be configured */
      if (nSM > EC_MAXSM)
         nSM = EC_MAXSM;
      /* iterate for every SM type defined */
      for (iSM = 2 ; iSM < nSM ; iSM++)
      {
         rdl = sizeof(tSM); tSM = 0;
         /* read SyncManager Communication Type */
         wkc = ecx_SDOread(context, Slave, ECT_SDO_SMCOMMTYPE, iSM + 1, FALSE, &rdl, &tSM, EC_TIMEOUTRXM);
         if (wkc > 0)
         {
// start slave bug prevention code, remove if possible
            if((iSM == 2) && (tSM == 2)) // SM2 has type 2 == mailbox out, this is a bug in the slave!
            {
               SMt_bug_add = 1; // try to correct, this works if the types are 0 1 2 3 and should be 1 2 3 4
            }
            if(tSM)
            {
               tSM += SMt_bug_add; // only add if SMt > 0
            }
            if((iSM == 2) && (tSM == 0)) // SM2 has type 0, this is a bug in the slave!
            {
               tSM = 3;
            }
            if((iSM == 3) && (tSM == 0)) // SM3 has type 0, this is a bug in the slave!
            {
               tSM = 4;
            }
// end slave bug prevention code

            context->slavelist[Slave].SMtype[iSM] = tSM;
            /* check if SM is unused -> clear enable flag */
            if (tSM == 0)
            {
               context->slavelist[Slave].SM[iSM].SMflags =
                  htoel( etohl(context->slavelist[Slave].SM[iSM].SMflags) & EC_SMENABLEMASK);
            }
            if ((tSM == 3) || (tSM == 4))
            {
               /* read the assign PDO */
               Tsize = ecx_readPDOassign(context, Slave, ECT_SDO_PDOASSIGN + iSM );
               /* if a mapping is found */
               if (Tsize)
               {
                  context->slavelist[Slave].SM[iSM].SMlength = htoes((Tsize + 7) / 8);
                  if (tSM == 3)
                  {
                     /* we are doing outputs */
                     *Osize += Tsize;
                  }
                  else
                  {
                     /* we are doing inputs */
                     *Isize += Tsize;
                  }
               }
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

/** CoE read PDO mapping in Complete Access mode (CA).
 *
 * CANopen has standard indexes defined for PDO mapping. This function
 * tries to read them and collect a full input and output mapping size
 * of designated slave. Slave has to support CA, otherwise use ec_readPDOmap().
 *
 * @param[in]  context  = context struct
 * @param[in]  Slave    = Slave number
 * @param[in]  Thread_n = Calling thread index
 * @param[out] Osize    = Size in bits of output mapping (rxPDO) found
 * @param[out] Isize    = Size in bits of input mapping (txPDO) found
 * @return >0 if mapping succesful.
 */
int ecx_readPDOmapCA(ecx_contextt *context, uint16 Slave, int Thread_n, int *Osize, int *Isize)
{
   int wkc, rdl;
   int retVal = 0;
   uint8 nSM, iSM, tSM;
   int Tsize;
   uint8 SMt_bug_add;

   *Isize = 0;
   *Osize = 0;
   SMt_bug_add = 0;
   rdl = sizeof(ec_SMcommtypet);
   context->SMcommtype[Thread_n].n = 0;
   /* read SyncManager Communication Type object count Complete Access*/
   wkc = ecx_SDOread(context, Slave, ECT_SDO_SMCOMMTYPE, 0x00, TRUE, &rdl,
         &(context->SMcommtype[Thread_n]), EC_TIMEOUTRXM);
   /* positive result from slave ? */
   if ((wkc > 0) && (context->SMcommtype[Thread_n].n > 2))
   {
      nSM = context->SMcommtype[Thread_n].n;
      /* limit to maximum number of SM defined, if true the slave can't be configured */
      if (nSM > EC_MAXSM)
      {
         nSM = EC_MAXSM;
         ecx_packeterror(context, Slave, 0, 0, 10); /* #SM larger than EC_MAXSM */
      }
      /* iterate for every SM type defined */
      for (iSM = 2 ; iSM < nSM ; iSM++)
      {
         tSM = context->SMcommtype[Thread_n].SMtype[iSM];

// start slave bug prevention code, remove if possible
         if((iSM == 2) && (tSM == 2)) // SM2 has type 2 == mailbox out, this is a bug in the slave!
         {
            SMt_bug_add = 1; // try to correct, this works if the types are 0 1 2 3 and should be 1 2 3 4
         }
         if(tSM)
         {
            tSM += SMt_bug_add; // only add if SMt > 0
         }
// end slave bug prevention code

         context->slavelist[Slave].SMtype[iSM] = tSM;
         /* check if SM is unused -> clear enable flag */
         if (tSM == 0)
         {
            context->slavelist[Slave].SM[iSM].SMflags =
               htoel( etohl(context->slavelist[Slave].SM[iSM].SMflags) & EC_SMENABLEMASK);
         }
         if ((tSM == 3) || (tSM == 4))
         {
            /* read the assign PDO */
            Tsize = ecx_readPDOassignCA(context, Slave, Thread_n,
                  ECT_SDO_PDOASSIGN + iSM );
            /* if a mapping is found */
            if (Tsize)
            {
               context->slavelist[Slave].SM[iSM].SMlength = htoes((Tsize + 7) / 8);
               if (tSM == 3)
               {
                  /* we are doing outputs */
                  *Osize += Tsize;
               }
               else
               {
                  /* we are doing inputs */
                  *Isize += Tsize;
               }
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

/** CoE read Object Description List.
 *
 * @param[in]  context  = context struct
 * @param[in]  Slave    = Slave number.
 * @param[out] pODlist  = resulting Object Description list.
 * @return Workcounter of slave response.
 */
int ecx_readODlist(ecx_contextt *context, uint16 Slave, ec_ODlistt *pODlist)
{
   ec_SDOservicet *SDOp, *aSDOp;
   ec_mbxbuft MbxIn, MbxOut;
   int wkc;
   uint16 x, n, i, sp, offset;
   boolean stop;
   uint8 cnt;
   boolean First;

   pODlist->Slave = Slave;
   pODlist->Entries = 0;
   ec_clearmbx(&MbxIn);
   /* clear pending out mailbox in slave if available. Timeout is set to 0 */
   wkc = ecx_mbxreceive(context, Slave, &MbxIn, 0);
   ec_clearmbx(&MbxOut);
   aSDOp = (ec_SDOservicet*)&MbxIn;
   SDOp = (ec_SDOservicet*)&MbxOut;
   SDOp->MbxHeader.length = htoes(0x0008);
   SDOp->MbxHeader.address = htoes(0x0000);
   SDOp->MbxHeader.priority = 0x00;
   /* Get new mailbox counter value */
   cnt = ec_nextmbxcnt(context->slavelist[Slave].mbx_cnt);
   context->slavelist[Slave].mbx_cnt = cnt;
   SDOp->MbxHeader.mbxtype = ECT_MBXT_COE + (cnt << 4); /* CoE */
   SDOp->CANOpen = htoes(0x000 + (ECT_COES_SDOINFO << 12)); /* number 9bits service upper 4 bits */
   SDOp->Opcode = ECT_GET_ODLIST_REQ; /* get object description list request */
   SDOp->Reserved = 0;
   SDOp->Fragments = 0; /* fragments left */
   SDOp->wdata[0] = htoes(0x01); /* all objects */
   /* send get object description list request to slave */
   wkc = ecx_mbxsend(context, Slave, &MbxOut, EC_TIMEOUTTXM);
   /* mailbox placed in slave ? */
   if (wkc > 0)
   {
      x = 0;
      sp = 0;
      First = TRUE;
      offset = 1; /* offset to skip info header in first frame, otherwise set to 0 */
      do
      {
         stop = TRUE; /* assume this is last iteration */
         ec_clearmbx(&MbxIn);
         /* read slave response */
         wkc = ecx_mbxreceive(context, Slave, &MbxIn, EC_TIMEOUTRXM);
         /* got response ? */
         if (wkc > 0)
         {
            /* response should be CoE and "get object description list response" */
            if (((aSDOp->MbxHeader.mbxtype & 0x0f) == ECT_MBXT_COE) &&
                ((aSDOp->Opcode & 0x7f) == ECT_GET_ODLIST_RES))
            {
               if (First)
               {
                  /* extract number of indexes from mailbox data size */
                  n = (etohs(aSDOp->MbxHeader.length) - (6 + 2)) / 2;
               }
               else
               {
                  /* extract number of indexes from mailbox data size */
                  n = (etohs(aSDOp->MbxHeader.length) - 6) / 2;
               }
               /* check if indexes fit in buffer structure */
               if ((sp + n) > EC_MAXODLIST)
               {
                  n = EC_MAXODLIST + 1 - sp;
                  ecx_SDOinfoerror(context, Slave, 0, 0, 0xf000000); /* Too many entries for master buffer */
                  stop = TRUE;
               }
               /* trim to maximum number of ODlist entries defined */
               if ((pODlist->Entries + n) > EC_MAXODLIST)
               {
                  n = EC_MAXODLIST - pODlist->Entries;
               }
               pODlist->Entries += n;
               /* extract indexes one by one */
               for (i = 0; i < n; i++)
               {
                  pODlist->Index[sp + i] = etohs(aSDOp->wdata[i + offset]);
               }
               sp += n;
               /* check if more fragments will follow */
               if (aSDOp->Fragments > 0)
               {
                  stop = FALSE;
               }
               First = FALSE;
               offset = 0;
            }
            /* got unexpected response from slave */
            else
            {
               if ((aSDOp->Opcode &  0x7f) == ECT_SDOINFO_ERROR) /* SDO info error received */
               {
                  ecx_SDOinfoerror(context, Slave, 0, 0, etohl(aSDOp->ldata[0]));
                  stop = TRUE;
               }
               else
               {
                  ecx_packeterror(context, Slave, 0, 0, 1); /* Unexpected frame returned */
               }
               wkc = 0;
               x += 20;
            }
         }
         x++;
      }
      while ((x <= 128) && !stop);
   }
   return wkc;
}

/** CoE read Object Description. Adds textual description to object indexes.
 *
 * @param[in]  context       = context struct
 * @param[in] Item           = Item number in ODlist.
 * @param[in,out] pODlist    = referencing Object Description list.
 * @return Workcounter of slave response.
 */
int ecx_readODdescription(ecx_contextt *context, uint16 Item, ec_ODlistt *pODlist)
{
   ec_SDOservicet *SDOp, *aSDOp;
   int wkc;
   uint16  n, Slave;
   ec_mbxbuft MbxIn, MbxOut;
   uint8 cnt;

   Slave = pODlist->Slave;
   pODlist->DataType[Item] = 0;
   pODlist->ObjectCode[Item] = 0;
   pODlist->MaxSub[Item] = 0;
   pODlist->Name[Item][0] = 0;
   ec_clearmbx(&MbxIn);
   /* clear pending out mailbox in slave if available. Timeout is set to 0 */
   wkc = ecx_mbxreceive(context, Slave, &MbxIn, 0);
   ec_clearmbx(&MbxOut);
   aSDOp = (ec_SDOservicet*)&MbxIn;
   SDOp = (ec_SDOservicet*)&MbxOut;
   SDOp->MbxHeader.length = htoes(0x0008);
   SDOp->MbxHeader.address = htoes(0x0000);
   SDOp->MbxHeader.priority = 0x00;
   /* Get new mailbox counter value */
   cnt = ec_nextmbxcnt(context->slavelist[Slave].mbx_cnt);
   context->slavelist[Slave].mbx_cnt = cnt;
   SDOp->MbxHeader.mbxtype = ECT_MBXT_COE + (cnt << 4); /* CoE */
   SDOp->CANOpen = htoes(0x000 + (ECT_COES_SDOINFO << 12)); /* number 9bits service upper 4 bits */
   SDOp->Opcode = ECT_GET_OD_REQ; /* get object description request */
   SDOp->Reserved = 0;
   SDOp->Fragments = 0; /* fragments left */
   SDOp->wdata[0] = htoes(pODlist->Index[Item]); /* Data of Index */
   /* send get object description request to slave */
   wkc = ecx_mbxsend(context, Slave, &MbxOut, EC_TIMEOUTTXM);
   /* mailbox placed in slave ? */
   if (wkc > 0)
   {
      ec_clearmbx(&MbxIn);
      /* read slave response */
      wkc = ecx_mbxreceive(context, Slave, &MbxIn, EC_TIMEOUTRXM);
      /* got response ? */
      if (wkc > 0)
      {
         if (((aSDOp->MbxHeader.mbxtype & 0x0f) == ECT_MBXT_COE) &&
             ((aSDOp->Opcode & 0x7f) == ECT_GET_OD_RES))
         {
            n = (etohs(aSDOp->MbxHeader.length) - 12); /* length of string(name of object) */
            if (n > EC_MAXNAME)
            {
               n = EC_MAXNAME; /* max chars */
            }
            pODlist->DataType[Item] = etohs(aSDOp->wdata[1]);
            pODlist->ObjectCode[Item] = aSDOp->bdata[5];
            pODlist->MaxSub[Item] = aSDOp->bdata[4];

            strncpy(pODlist->Name[Item] , (char *)&aSDOp->bdata[6], n);
            pODlist->Name[Item][n] = 0x00; /* String terminator */
         }
         /* got unexpected response from slave */
         else
         {
            if (((aSDOp->Opcode & 0x7f) == ECT_SDOINFO_ERROR)) /* SDO info error received */
            {
               ecx_SDOinfoerror(context, Slave,pODlist->Index[Item], 0, etohl(aSDOp->ldata[0]));
            }
            else
            {
               ecx_packeterror(context, Slave,pODlist->Index[Item], 0, 1); /* Unexpected frame returned */
            }
            wkc = 0;
         }
      }
   }

   return wkc;
}

/** CoE read SDO service object entry, single subindex.
 * Used in ec_readOE().
 *
 * @param[in]  context       = context struct
 * @param[in] Item           = Item in ODlist.
 * @param[in] SubI           = Subindex of item in ODlist.
 * @param[in] pODlist        = Object description list for reference.
 * @param[out] pOElist       = resulting object entry structure.
 * @return Workcounter of slave response.
 */
int ecx_readOEsingle(ecx_contextt *context, uint16 Item, uint8 SubI, ec_ODlistt *pODlist, ec_OElistt *pOElist)
{
   ec_SDOservicet *SDOp, *aSDOp;
   int wkc;
   uint16 Index, Slave;
   int16 n;
   ec_mbxbuft MbxIn, MbxOut;
   uint8 cnt;

   wkc = 0;
   Slave = pODlist->Slave;
   Index = pODlist->Index[Item];
   ec_clearmbx(&MbxIn);
   /* clear pending out mailbox in slave if available. Timeout is set to 0 */
   wkc = ecx_mbxreceive(context, Slave, &MbxIn, 0);
   ec_clearmbx(&MbxOut);
   aSDOp = (ec_SDOservicet*)&MbxIn;
   SDOp = (ec_SDOservicet*)&MbxOut;
   SDOp->MbxHeader.length = htoes(0x000a);
   SDOp->MbxHeader.address = htoes(0x0000);
   SDOp->MbxHeader.priority = 0x00;
   /* Get new mailbox counter value */
   cnt = ec_nextmbxcnt(context->slavelist[Slave].mbx_cnt);
   context->slavelist[Slave].mbx_cnt = cnt;
   SDOp->MbxHeader.mbxtype = ECT_MBXT_COE + (cnt << 4); /* CoE */
   SDOp->CANOpen = htoes(0x000 + (ECT_COES_SDOINFO << 12)); /* number 9bits service upper 4 bits */
   SDOp->Opcode = ECT_GET_OE_REQ; /* get object entry description request */
   SDOp->Reserved = 0;
   SDOp->Fragments = 0;      /* fragments left */
   SDOp->wdata[0] = htoes(Index);      /* Index */
   SDOp->bdata[2] = SubI;       /* SubIndex */
   SDOp->bdata[3] = 1 + 2 + 4; /* get access rights, object category, PDO */
   /* send get object entry description request to slave */
   wkc = ecx_mbxsend(context, Slave, &MbxOut, EC_TIMEOUTTXM);
   /* mailbox placed in slave ? */
   if (wkc > 0)
   {
      ec_clearmbx(&MbxIn);
      /* read slave response */
      wkc = ecx_mbxreceive(context, Slave, &MbxIn, EC_TIMEOUTRXM);
      /* got response ? */
      if (wkc > 0)
      {
         if (((aSDOp->MbxHeader.mbxtype & 0x0f) == ECT_MBXT_COE) &&
             ((aSDOp->Opcode &  0x7f) == ECT_GET_OE_RES))
         {
            pOElist->Entries++;
            n = (etohs(aSDOp->MbxHeader.length) - 16); /* length of string(name of object) */
            if (n > EC_MAXNAME)
            {
               n = EC_MAXNAME; /* max string length */
            }
            if (n < 0 )
            {
               n = 0;
            }
            pOElist->ValueInfo[SubI] = aSDOp->bdata[3];
            pOElist->DataType[SubI] = etohs(aSDOp->wdata[2]);
            pOElist->BitLength[SubI] = etohs(aSDOp->wdata[3]);
            pOElist->ObjAccess[SubI] = etohs(aSDOp->wdata[4]);

            strncpy(pOElist->Name[SubI] , (char *)&aSDOp->wdata[5], n);
            pOElist->Name[SubI][n] = 0x00; /* string terminator */
         }
         /* got unexpected response from slave */
         else
         {
            if (((aSDOp->Opcode & 0x7f) == ECT_SDOINFO_ERROR)) /* SDO info error received */
            {
               ecx_SDOinfoerror(context, Slave, Index, SubI, etohl(aSDOp->ldata[0]));
            }
            else
            {
               ecx_packeterror(context, Slave, Index, SubI, 1); /* Unexpected frame returned */
            }
            wkc = 0;
         }
      }
   }

   return wkc;
}

/** CoE read SDO service object entry.
 *
 * @param[in] context        = context struct
 * @param[in] Item           = Item in ODlist.
 * @param[in] pODlist        = Object description list for reference.
 * @param[out] pOElist       = resulting object entry structure.
 * @return Workcounter of slave response.
 */
int ecx_readOE(ecx_contextt *context, uint16 Item, ec_ODlistt *pODlist, ec_OElistt *pOElist)
{
   uint16 SubCount;
   int wkc;
   uint8 SubI;

   wkc = 0;
   pOElist->Entries = 0;
   SubI = pODlist->MaxSub[Item];
   /* for each entry found in ODlist */
   for (SubCount = 0; SubCount <= SubI; SubCount++)
   {
      /* read subindex of entry */
      wkc = ecx_readOEsingle(context, Item, (uint8)SubCount, pODlist, pOElist);
   }

   return wkc;
}

#ifdef EC_VER1
/** Report SDO error.
 *
 * @param[in]  Slave      = Slave number
 * @param[in]  Index      = Index that generated error
 * @param[in]  SubIdx     = Subindex that generated error
 * @param[in]  AbortCode  = Abortcode, see EtherCAT documentation for list
 * @see ecx_SDOerror
 */
void ec_SDOerror(uint16 Slave, uint16 Index, uint8 SubIdx, int32 AbortCode)
{
   ecx_SDOerror(&ecx_context, Slave, Index, SubIdx, AbortCode);
}

/** CoE SDO read, blocking. Single subindex or Complete Access.
 *
 * Only a "normal" upload request is issued. If the requested parameter is <= 4bytes
 * then a "expedited" response is returned, otherwise a "normal" response. If a "normal"
 * response is larger than the mailbox size then the response is segmented. The function
 * will combine all segments and copy them to the parameter buffer.
 *
 * @param[in]  slave      = Slave number
 * @param[in]  index      = Index to read
 * @param[in]  subindex   = Subindex to read, must be 0 or 1 if CA is used.
 * @param[in]  CA         = FALSE = single subindex. TRUE = Complete Access, all subindexes read.
 * @param[in,out] psize   = Size in bytes of parameter buffer, returns bytes read from SDO.
 * @param[out] p          = Pointer to parameter buffer
 * @param[in]  timeout    = Timeout in us, standard is EC_TIMEOUTRXM
 * @return Workcounter from last slave response
 * @see ecx_SDOread
 */
int ec_SDOread(uint16 slave, uint16 index, uint8 subindex,
               boolean CA, int *psize, void *p, int timeout)
{
   return ecx_SDOread(&ecx_context, slave, index, subindex, CA, psize, p, timeout);
}

/** CoE SDO write, blocking. Single subindex or Complete Access.
 *
 * A "normal" download request is issued, unless we have
 * small data, then a "expedited" transfer is used. If the parameter is larger than
 * the mailbox size then the download is segmented. The function will split the
 * parameter data in segments and send them to the slave one by one.
 *
 * @param[in]  Slave      = Slave number
 * @param[in]  Index      = Index to write
 * @param[in]  SubIndex   = Subindex to write, must be 0 or 1 if CA is used.
 * @param[in]  CA         = FALSE = single subindex. TRUE = Complete Access, all subindexes written.
 * @param[in]  psize      = Size in bytes of parameter buffer.
 * @param[out] p          = Pointer to parameter buffer
 * @param[in]  Timeout    = Timeout in us, standard is EC_TIMEOUTRXM
 * @return Workcounter from last slave response
 * @see ecx_SDOwrite
 */
int ec_SDOwrite(uint16 Slave, uint16 Index, uint8 SubIndex,
                boolean CA, int psize, void *p, int Timeout)
{
   return ecx_SDOwrite(&ecx_context, Slave, Index, SubIndex, CA, psize, p, Timeout);
}

/** CoE RxPDO write, blocking.
 *
 * A RxPDO download request is issued.
 *
 * @param[in]  Slave         = Slave number
 * @param[in]  RxPDOnumber   = Related RxPDO number
 * @param[in]  psize         = Size in bytes of PDO buffer.
 * @param[out] p             = Pointer to PDO buffer
 * @return Workcounter from last slave response
 * @see ecx_RxPDO
 */
int ec_RxPDO(uint16 Slave, uint16 RxPDOnumber, int psize, void *p)
{
   return ecx_RxPDO(&ecx_context, Slave, RxPDOnumber, psize, p);
}

/** CoE TxPDO read remote request, blocking.
 *
 * A RxPDO download request is issued.
 *
 * @param[in]  slave         = Slave number
 * @param[in]  TxPDOnumber   = Related TxPDO number
 * @param[in,out] psize      = Size in bytes of PDO buffer, returns bytes read from PDO.
 * @param[out] p             = Pointer to PDO buffer
 * @param[in]  timeout       = Timeout in us, standard is EC_TIMEOUTRXM
 * @return Workcounter from last slave response
 * @see ecx_TxPDO
 */
int ec_TxPDO(uint16 slave, uint16 TxPDOnumber , int *psize, void *p, int timeout)
{
   return ecx_TxPDO(&ecx_context, slave, TxPDOnumber, psize, p, timeout);
}

/** Read PDO assign structure
 * @param[in]  Slave         = Slave number
 * @param[in]  PDOassign     = PDO assign object
 * @return total bitlength of PDO assign
 */
int ec_readPDOassign(uint16 Slave, uint16 PDOassign)
{
   return ecx_readPDOassign(&ecx_context, Slave, PDOassign);
}

/** Read PDO assign structure in Complete Access mode
 * @param[in]  Slave         = Slave number
 * @param[in]  PDOassign     = PDO assign object
 * @param[in]  Thread_n      = Calling thread index
 * @return total bitlength of PDO assign
 * @see ecx_readPDOmap
 */
int ec_readPDOassignCA(uint16 Slave, uint16 PDOassign, int Thread_n)
{
   return ecx_readPDOassignCA(&ecx_context, Slave, Thread_n, PDOassign);
}

/** CoE read PDO mapping.
 *
 * CANopen has standard indexes defined for PDO mapping. This function
 * tries to read them and collect a full input and output mapping size
 * of designated slave.
 *
 * For details, see #ecx_readPDOmap
 *
 * @param[in] Slave    = Slave number
 * @param[out] Osize   = Size in bits of output mapping (rxPDO) found
 * @param[out] Isize   = Size in bits of input mapping (txPDO) found
 * @return >0 if mapping succesful.
 */
int ec_readPDOmap(uint16 Slave, int *Osize, int *Isize)
{
   return ecx_readPDOmap(&ecx_context, Slave, Osize, Isize);
}

/** CoE read PDO mapping in Complete Access mode (CA).
 *
 * CANopen has standard indexes defined for PDO mapping. This function
 * tries to read them and collect a full input and output mapping size
 * of designated slave. Slave has to support CA, otherwise use ec_readPDOmap().
 *
 * @param[in] Slave    = Slave number
 * @param[in] Thread_n = Calling thread index
 * @param[out] Osize   = Size in bits of output mapping (rxPDO) found
 * @param[out] Isize   = Size in bits of input mapping (txPDO) found
 * @return >0 if mapping succesful.
 * @see ecx_readPDOmap ec_readPDOmapCA
 */
int ec_readPDOmapCA(uint16 Slave, int Thread_n, int *Osize, int *Isize)
{
   return ecx_readPDOmapCA(&ecx_context, Slave, Thread_n, Osize, Isize);
}

/** CoE read Object Description List.
 *
 * @param[in] Slave      = Slave number.
 * @param[out] pODlist  = resulting Object Description list.
 * @return Workcounter of slave response.
 * @see ecx_readODlist
 */
int ec_readODlist(uint16 Slave, ec_ODlistt *pODlist)
{
   return ecx_readODlist(&ecx_context, Slave, pODlist);
}

/** CoE read Object Description. Adds textual description to object indexes.
 *
 * @param[in] Item           = Item number in ODlist.
 * @param[in,out] pODlist    = referencing Object Description list.
 * @return Workcounter of slave response.
 * @see ecx_readODdescription
 */
int ec_readODdescription(uint16 Item, ec_ODlistt *pODlist)
{
   return ecx_readODdescription(&ecx_context, Item, pODlist);
}

int ec_readOEsingle(uint16 Item, uint8 SubI, ec_ODlistt *pODlist, ec_OElistt *pOElist)
{
   return ecx_readOEsingle(&ecx_context, Item, SubI, pODlist, pOElist);
}

/** CoE read SDO service object entry.
 *
 * @param[in] Item           = Item in ODlist.
 * @param[in] pODlist        = Object description list for reference.
 * @param[out] pOElist       = resulting object entry structure.
 * @return Workcounter of slave response.
 * @see ecx_readOE
 */
int ec_readOE(uint16 Item, ec_ODlistt *pODlist, ec_OElistt *pOElist)
{
   return ecx_readOE(&ecx_context, Item, pODlist, pOElist);
}
#endif
