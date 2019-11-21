/*
 * Simple Open EtherCAT Master Library
 *
 * File    : ethercatmain.c
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

/**
 * \file
 * \brief
 * Main EtherCAT functions.
 *
 * Initialisation, state set and read, mailbox primitives, EEPROM primitives,
 * SII reading and processdata exchange.
 *
 * Defines ec_slave[]. All slave information is put in this structure.
 * Needed for most user interaction with slaves.
 */

#include <stdio.h>
#include <string.h>
#include "osal.h"
#include "oshw.h"
#include "ethercattype.h"
#include "ethercatbase.h"
#include "ethercatmain.h"


/** delay in us for eeprom ready loop */
#define EC_LOCALDELAY  200

/** record for ethercat eeprom communications */
PACKED_BEGIN
typedef struct PACKED
{
   uint16    comm;
   uint16    addr;
   uint16    d2;
} ec_eepromt;
PACKED_END

/** mailbox error structure */
PACKED_BEGIN
typedef struct PACKED
{
   ec_mbxheadert   MbxHeader;
   uint16          Type;
   uint16          Detail;
} ec_mbxerrort;
PACKED_END

/** emergency request structure */
PACKED_BEGIN
typedef struct PACKED
{
   ec_mbxheadert   MbxHeader;
   uint16          CANOpen;
   uint16          ErrorCode;
   uint8           ErrorReg;
   uint8           bData;
   uint16          w1,w2;
} ec_emcyt;
PACKED_END

#ifdef EC_VER1
/** Main slave data array.
 *  Each slave found on the network gets its own record.
 *  ec_slave[0] is reserved for the master. Structure gets filled
 *  in by the configuration function ec_config().
 */
ec_slavet               ec_slave[EC_MAXSLAVE];
/** number of slaves found on the network */
int                     ec_slavecount;
/** slave group structure */
ec_groupt               ec_group[EC_MAXGROUP];

/** cache for EEPROM read functions */
static uint8            ec_esibuf[EC_MAXEEPBUF];
/** bitmap for filled cache buffer bytes */
static uint32           ec_esimap[EC_MAXEEPBITMAP];
/** current slave for EEPROM cache buffer */
static ec_eringt        ec_elist;
static ec_idxstackT     ec_idxstack;

/** SyncManager Communication Type struct to store data of one slave */
static ec_SMcommtypet   ec_SMcommtype[EC_MAX_MAPT];
/** PDO assign struct to store data of one slave */
static ec_PDOassignt    ec_PDOassign[EC_MAX_MAPT];
/** PDO description struct to store data of one slave */
static ec_PDOdesct      ec_PDOdesc[EC_MAX_MAPT];

/** buffer for EEPROM SM data */
static ec_eepromSMt     ec_SM;
/** buffer for EEPROM FMMU data */
static ec_eepromFMMUt   ec_FMMU;
/** Global variable TRUE if error available in error stack */
boolean                 EcatError = FALSE;

int64                   ec_DCtime;

ecx_portt               ecx_port;
ecx_redportt            ecx_redport;

ecx_contextt  ecx_context = {
    &ecx_port,          // .port          =
    &ec_slave[0],       // .slavelist     =
    &ec_slavecount,     // .slavecount    =
    EC_MAXSLAVE,        // .maxslave      =
    &ec_group[0],       // .grouplist     =
    EC_MAXGROUP,        // .maxgroup      =
    &ec_esibuf[0],      // .esibuf        =
    &ec_esimap[0],      // .esimap        =
    0,                  // .esislave      =
    &ec_elist,          // .elist         =
    &ec_idxstack,       // .idxstack      =
    &EcatError,         // .ecaterror     =
    0,                  // .DCtO          =
    0,                  // .DCl           =
    &ec_DCtime,         // .DCtime        =
    &ec_SMcommtype[0],  // .SMcommtype    =
    &ec_PDOassign[0],   // .PDOassign     =
    &ec_PDOdesc[0],     // .PDOdesc       =
    &ec_SM,             // .eepSM         =
    &ec_FMMU,           // .eepFMMU       =
    NULL                // .FOEhook()
};
#endif

/** Create list over available network adapters.
 *
 * @return First element in list over available network adapters.
 */
ec_adaptert * ec_find_adapters (void)
{
   ec_adaptert * ret_adapter;

   ret_adapter = oshw_find_adapters ();

   return ret_adapter;
}

/** Free dynamically allocated list over available network adapters.
 *
 * @param[in] adapter = Struct holding adapter name, description and pointer to next.
 */
void ec_free_adapters (ec_adaptert * adapter)
{
   oshw_free_adapters (adapter);
}

/** Pushes an error on the error list.
 *
 * @param[in] context        = context struct
 * @param[in] Ec pointer describing the error.
 */
void ecx_pusherror(ecx_contextt *context, const ec_errort *Ec)
{
   context->elist->Error[context->elist->head] = *Ec;
   context->elist->Error[context->elist->head].Signal = TRUE;
   context->elist->head++;
   if (context->elist->head > EC_MAXELIST)
   {
      context->elist->head = 0;
   }
   if (context->elist->head == context->elist->tail)
   {
      context->elist->tail++;
   }
   if (context->elist->tail > EC_MAXELIST)
   {
      context->elist->tail = 0;
   }
   *(context->ecaterror) = TRUE;
}

/** Pops an error from the list.
 *
 * @param[in] context        = context struct
 * @param[out] Ec = Struct describing the error.
 * @return TRUE if an error was popped.
 */
boolean ecx_poperror(ecx_contextt *context, ec_errort *Ec)
{
   boolean notEmpty = (context->elist->head != context->elist->tail);

   *Ec = context->elist->Error[context->elist->tail];
   context->elist->Error[context->elist->tail].Signal = FALSE;
   if (notEmpty)
   {
      context->elist->tail++;
      if (context->elist->tail > EC_MAXELIST)
      {
         context->elist->tail = 0;
      }
   }
   else
   {
      *(context->ecaterror) = FALSE;
   }
   return notEmpty;
}

/** Check if error list has entries.
 *
 * @param[in] context        = context struct
 * @return TRUE if error list contains entries.
 */
boolean ecx_iserror(ecx_contextt *context)
{
   return (context->elist->head != context->elist->tail);
}

/** Report packet error
 *
 * @param[in]  context        = context struct
 * @param[in]  Slave      = Slave number
 * @param[in]  Index      = Index that generated error
 * @param[in]  SubIdx     = Subindex that generated error
 * @param[in]  ErrorCode  = Error code
 */
void ecx_packeterror(ecx_contextt *context, uint16 Slave, uint16 Index, uint8 SubIdx, uint16 ErrorCode)
{
   ec_errort Ec;

   memset(&Ec, 0, sizeof(Ec));
   Ec.Time = osal_current_time();
   Ec.Slave = Slave;
   Ec.Index = Index;
   Ec.SubIdx = SubIdx;
   *(context->ecaterror) = TRUE;
   Ec.Etype = EC_ERR_TYPE_PACKET_ERROR;
   Ec.ErrorCode = ErrorCode;
   ecx_pusherror(context, &Ec);
}

/** Report Mailbox Error
 *
 * @param[in]  context        = context struct
 * @param[in]  Slave        = Slave number
 * @param[in]  Detail       = Following EtherCAT specification
 */
static void ecx_mbxerror(ecx_contextt *context, uint16 Slave,uint16 Detail)
{
   ec_errort Ec;

   memset(&Ec, 0, sizeof(Ec));
   Ec.Time = osal_current_time();
   Ec.Slave = Slave;
   Ec.Index = 0;
   Ec.SubIdx = 0;
   Ec.Etype = EC_ERR_TYPE_MBX_ERROR;
   Ec.ErrorCode = Detail;
   ecx_pusherror(context, &Ec);
}

/** Report Mailbox Emergency Error
 *
 * @param[in]  context        = context struct
 * @param[in]  Slave      = Slave number
 * @param[in]  ErrorCode  = Following EtherCAT specification
 * @param[in]  ErrorReg
 * @param[in]  b1
 * @param[in]  w1
 * @param[in]  w2
 */
static void ecx_mbxemergencyerror(ecx_contextt *context, uint16 Slave,uint16 ErrorCode,uint16 ErrorReg,
    uint8 b1, uint16 w1, uint16 w2)
{
   ec_errort Ec;

   memset(&Ec, 0, sizeof(Ec));
   Ec.Time = osal_current_time();
   Ec.Slave = Slave;
   Ec.Index = 0;
   Ec.SubIdx = 0;
   Ec.Etype = EC_ERR_TYPE_EMERGENCY;
   Ec.ErrorCode = ErrorCode;
   Ec.ErrorReg = (uint8)ErrorReg;
   Ec.b1 = b1;
   Ec.w1 = w1;
   Ec.w2 = w2;
   ecx_pusherror(context, &Ec);
}

/** Initialise lib in single NIC mode
 * @param[in]  context = context struct
 * @param[in] ifname   = Dev name, f.e. "eth0"
 * @return >0 if OK
 */
int ecx_init(ecx_contextt *context, const char * ifname)
{
   return ecx_setupnic(context->port, ifname, FALSE);
}

/** Initialise lib in redundant NIC mode
 * @param[in]  context  = context struct
 * @param[in]  redport  = pointer to redport, redundant port data
 * @param[in]  ifname   = Primary Dev name, f.e. "eth0"
 * @param[in]  if2name  = Secondary Dev name, f.e. "eth1"
 * @return >0 if OK
 */
int ecx_init_redundant(ecx_contextt *context, ecx_redportt *redport, const char *ifname, char *if2name)
{
   int rval, zbuf;
   ec_etherheadert *ehp;

   context->port->redport = redport;
   ecx_setupnic(context->port, ifname, FALSE);
   rval = ecx_setupnic(context->port, if2name, TRUE);
   /* prepare "dummy" BRD tx frame for redundant operation */
   ehp = (ec_etherheadert *)&(context->port->txbuf2);
   ehp->sa1 = oshw_htons(secMAC[0]);
   zbuf = 0;
   ecx_setupdatagram(context->port, &(context->port->txbuf2), EC_CMD_BRD, 0, 0x0000, 0x0000, 2, &zbuf);
   context->port->txbuflength2 = ETH_HEADERSIZE + EC_HEADERSIZE + EC_WKCSIZE + 2;

   return rval;
}

/** Close lib.
 * @param[in]  context        = context struct
 */
void ecx_close(ecx_contextt *context)
{
   ecx_closenic(context->port);
};

/** Read one byte from slave EEPROM via cache.
 *  If the cache location is empty then a read request is made to the slave.
 *  Depending on the slave capabillities the request is 4 or 8 bytes.
 *  @param[in] context = context struct
 *  @param[in] slave   = slave number
 *  @param[in] address = eeprom address in bytes (slave uses words)
 *  @return requested byte, if not available then 0xff
 */
uint8 ecx_siigetbyte(ecx_contextt *context, uint16 slave, uint16 address)
{
   uint16 configadr, eadr;
   uint64 edat;
   uint16 mapw, mapb;
   int lp,cnt;
   uint8 retval;

   retval = 0xff;
   if (slave != context->esislave) /* not the same slave? */
   {
      memset(context->esimap, 0x00, EC_MAXEEPBITMAP * sizeof(uint32)); /* clear esibuf cache map */
      context->esislave = slave;
   }
   if (address < EC_MAXEEPBUF)
   {
      mapw = address >> 5;
      mapb = address - (mapw << 5);
      if (context->esimap[mapw] & (uint32)(1 << mapb))
      {
         /* byte is already in buffer */
         retval = context->esibuf[address];
      }
      else
      {
         /* byte is not in buffer, put it there */
         configadr = context->slavelist[slave].configadr;
         ecx_eeprom2master(context, slave); /* set eeprom control to master */
         eadr = address >> 1;
         edat = ecx_readeepromFP (context, configadr, eadr, EC_TIMEOUTEEP);
         /* 8 byte response */
         if (context->slavelist[slave].eep_8byte)
         {
            put_unaligned64(edat, &(context->esibuf[eadr << 1]));
            cnt = 8;
         }
         /* 4 byte response */
         else
         {
            put_unaligned32(edat, &(context->esibuf[eadr << 1]));
            cnt = 4;
         }
         /* find bitmap location */
         mapw = eadr >> 4;
         mapb = (eadr << 1) - (mapw << 5);
         for(lp = 0 ; lp < cnt ; lp++)
         {
            /* set bitmap for each byte that is read */
            context->esimap[mapw] |= (1 << mapb);
            mapb++;
            if (mapb > 31)
            {
               mapb = 0;
               mapw++;
            }
         }
         retval = context->esibuf[address];
      }
   }

   return retval;
}

/** Find SII section header in slave EEPROM.
 *  @param[in]  context        = context struct
 *  @param[in] slave   = slave number
 *  @param[in] cat     = section category
 *  @return byte address of section at section length entry, if not available then 0
 */
int16 ecx_siifind(ecx_contextt *context, uint16 slave, uint16 cat)
{
   int16 a;
   uint16 p;
   uint8 eectl = context->slavelist[slave].eep_pdi;

   a = ECT_SII_START << 1;
   /* read first SII section category */
   p = ecx_siigetbyte(context, slave, a++);
   p += (ecx_siigetbyte(context, slave, a++) << 8);
   /* traverse SII while category is not found and not EOF */
   while ((p != cat) && (p != 0xffff))
   {
      /* read section length */
      p = ecx_siigetbyte(context, slave, a++);
      p += (ecx_siigetbyte(context, slave, a++) << 8);
      /* locate next section category */
      a += p << 1;
      /* read section category */
      p = ecx_siigetbyte(context, slave, a++);
      p += (ecx_siigetbyte(context, slave, a++) << 8);
   }
   if (p != cat)
   {
      a = 0;
   }
   if (eectl)
   {
      ecx_eeprom2pdi(context, slave); /* if eeprom control was previously pdi then restore */
   }

   return a;
}

/** Get string from SII string section in slave EEPROM.
 *  @param[in]  context = context struct
 *  @param[out] str     = requested string, 0x00 if not found
 *  @param[in]  slave   = slave number
 *  @param[in]  Sn      = string number
 */
void ecx_siistring(ecx_contextt *context, char *str, uint16 slave, uint16 Sn)
{
   uint16 a,i,j,l,n,ba;
   char *ptr;
   uint8 eectl = context->slavelist[slave].eep_pdi;

   ptr = str;
   a = ecx_siifind (context, slave, ECT_SII_STRING); /* find string section */
   if (a > 0)
   {
      ba = a + 2; /* skip SII section header */
      n = ecx_siigetbyte(context, slave, ba++); /* read number of strings in section */
      if (Sn <= n) /* is req string available? */
      {
         for (i = 1; i <= Sn; i++) /* walk through strings */
         {
            l = ecx_siigetbyte(context, slave, ba++); /* length of this string */
            if (i < Sn)
            {
               ba += l;
            }
            else
            {
               ptr = str;
               for (j = 1; j <= l; j++) /* copy one string */
               {
                  if(j <= EC_MAXNAME)
                  {
                     *ptr = (char)ecx_siigetbyte(context, slave, ba++);
                     ptr++;
                  }
                  else
                  {
                     ba++;
                  }
               }
            }
         }
         *ptr = 0; /* add zero terminator */
      }
      else
      {
         ptr = str;
         *ptr = 0; /* empty string */
      }
   }
   if (eectl)
   {
      ecx_eeprom2pdi(context, slave); /* if eeprom control was previously pdi then restore */
   }
}

/** Get FMMU data from SII FMMU section in slave EEPROM.
 *  @param[in]  context = context struct
 *  @param[in]  slave   = slave number
 *  @param[out] FMMU    = FMMU struct from SII, max. 4 FMMU's
 *  @return number of FMMU's defined in section
 */
uint16 ecx_siiFMMU(ecx_contextt *context, uint16 slave, ec_eepromFMMUt* FMMU)
{
   uint16  a;
   uint8 eectl = context->slavelist[slave].eep_pdi;

   FMMU->nFMMU = 0;
   FMMU->FMMU0 = 0;
   FMMU->FMMU1 = 0;
   FMMU->FMMU2 = 0;
   FMMU->FMMU3 = 0;
   FMMU->Startpos = ecx_siifind(context, slave, ECT_SII_FMMU);

   if (FMMU->Startpos > 0)
   {
      a = FMMU->Startpos;
      FMMU->nFMMU = ecx_siigetbyte(context, slave, a++);
      FMMU->nFMMU += (ecx_siigetbyte(context, slave, a++) << 8);
      FMMU->nFMMU *= 2;
      FMMU->FMMU0 = ecx_siigetbyte(context, slave, a++);
      FMMU->FMMU1 = ecx_siigetbyte(context, slave, a++);
      if (FMMU->nFMMU > 2)
      {
         FMMU->FMMU2 = ecx_siigetbyte(context, slave, a++);
         FMMU->FMMU3 = ecx_siigetbyte(context, slave, a++);
      }
   }
   if (eectl)
   {
      ecx_eeprom2pdi(context, slave); /* if eeprom control was previously pdi then restore */
   }

   return FMMU->nFMMU;
}

/** Get SM data from SII SM section in slave EEPROM.
 *  @param[in]  context = context struct
 *  @param[in]  slave   = slave number
 *  @param[out] SM      = first SM struct from SII
 *  @return number of SM's defined in section
 */
uint16 ecx_siiSM(ecx_contextt *context, uint16 slave, ec_eepromSMt* SM)
{
   uint16 a,w;
   uint8 eectl = context->slavelist[slave].eep_pdi;

   SM->nSM = 0;
   SM->Startpos = ecx_siifind(context, slave, ECT_SII_SM);
   if (SM->Startpos > 0)
   {
      a = SM->Startpos;
      w = ecx_siigetbyte(context, slave, a++);
      w += (ecx_siigetbyte(context, slave, a++) << 8);
      SM->nSM = (w / 4);
      SM->PhStart = ecx_siigetbyte(context, slave, a++);
      SM->PhStart += (ecx_siigetbyte(context, slave, a++) << 8);
      SM->Plength = ecx_siigetbyte(context, slave, a++);
      SM->Plength += (ecx_siigetbyte(context, slave, a++) << 8);
      SM->Creg = ecx_siigetbyte(context, slave, a++);
      SM->Sreg = ecx_siigetbyte(context, slave, a++);
      SM->Activate = ecx_siigetbyte(context, slave, a++);
      SM->PDIctrl = ecx_siigetbyte(context, slave, a++);
   }
   if (eectl)
   {
      ecx_eeprom2pdi(context, slave); /* if eeprom control was previously pdi then restore */
   }

   return SM->nSM;
}

/** Get next SM data from SII SM section in slave EEPROM.
 *  @param[in]  context = context struct
 *  @param[in]  slave   = slave number
 *  @param[out] SM      = first SM struct from SII
 *  @param[in]  n       = SM number
 *  @return >0 if OK
 */
uint16 ecx_siiSMnext(ecx_contextt *context, uint16 slave, ec_eepromSMt* SM, uint16 n)
{
   uint16 a;
   uint16 retVal = 0;
   uint8 eectl = context->slavelist[slave].eep_pdi;

   if (n < SM->nSM)
   {
      a = SM->Startpos + 2 + (n * 8);
      SM->PhStart = ecx_siigetbyte(context, slave, a++);
      SM->PhStart += (ecx_siigetbyte(context, slave, a++) << 8);
      SM->Plength = ecx_siigetbyte(context, slave, a++);
      SM->Plength += (ecx_siigetbyte(context, slave, a++) << 8);
      SM->Creg = ecx_siigetbyte(context, slave, a++);
      SM->Sreg = ecx_siigetbyte(context, slave, a++);
      SM->Activate = ecx_siigetbyte(context, slave, a++);
      SM->PDIctrl = ecx_siigetbyte(context, slave, a++);
      retVal = 1;
   }
   if (eectl)
   {
      ecx_eeprom2pdi(context, slave); /* if eeprom control was previously pdi then restore */
   }

   return retVal;
}

/** Get PDO data from SII PDO section in slave EEPROM.
 *  @param[in]  context = context struct
 *  @param[in]  slave   = slave number
 *  @param[out] PDO     = PDO struct from SII
 *  @param[in]  t       = 0=RXPDO 1=TXPDO
 *  @return mapping size in bits of PDO
 */
int ecx_siiPDO(ecx_contextt *context, uint16 slave, ec_eepromPDOt* PDO, uint8 t)
{
   uint16 a , w, c, e, er, Size;
   uint8 eectl = context->slavelist[slave].eep_pdi;

   Size = 0;
   PDO->nPDO = 0;
   PDO->Length = 0;
   PDO->Index[1] = 0;
   for (c = 0 ; c < EC_MAXSM ; c++) PDO->SMbitsize[c] = 0;
   if (t > 1)
      t = 1;
   PDO->Startpos = ecx_siifind(context, slave, ECT_SII_PDO + t);
   if (PDO->Startpos > 0)
   {
      a = PDO->Startpos;
      w = ecx_siigetbyte(context, slave, a++);
      w += (ecx_siigetbyte(context, slave, a++) << 8);
      PDO->Length = w;
      c = 1;
      /* traverse through all PDOs */
      do
      {
         PDO->nPDO++;
         PDO->Index[PDO->nPDO] = ecx_siigetbyte(context, slave, a++);
         PDO->Index[PDO->nPDO] += (ecx_siigetbyte(context, slave, a++) << 8);
         PDO->BitSize[PDO->nPDO] = 0;
         c++;
         e = ecx_siigetbyte(context, slave, a++);
         PDO->SyncM[PDO->nPDO] = ecx_siigetbyte(context, slave, a++);
         a += 4;
         c += 2;
         if (PDO->SyncM[PDO->nPDO] < EC_MAXSM) /* active and in range SM? */
         {
            /* read all entries defined in PDO */
            for (er = 1; er <= e; er++)
            {
               c += 4;
               a += 5;
               PDO->BitSize[PDO->nPDO] += ecx_siigetbyte(context, slave, a++);
               a += 2;
            }
            PDO->SMbitsize[ PDO->SyncM[PDO->nPDO] ] += PDO->BitSize[PDO->nPDO];
            Size += PDO->BitSize[PDO->nPDO];
            c++;
         }
         else /* PDO deactivated because SM is 0xff or > EC_MAXSM */
         {
            c += 4 * e;
            a += 8 * e;
            c++;
         }
         if (PDO->nPDO >= (EC_MAXEEPDO - 1))
         {
            c = PDO->Length; /* limit number of PDO entries in buffer */
         }
      }
      while (c < PDO->Length);
   }
   if (eectl)
   {
      ecx_eeprom2pdi(context, slave); /* if eeprom control was previously pdi then restore */
   }

   return (Size);
}

#define MAX_FPRD_MULTI 64

int ecx_FPRD_multi(ecx_contextt *context, int n, uint16 *configlst, ec_alstatust *slstatlst, int timeout)
{
   int wkc;
   uint8 idx;
   ecx_portt *port;
   int sldatapos[MAX_FPRD_MULTI];
   int slcnt;

   port = context->port;
   idx = ecx_getindex(port);
   slcnt = 0;
   ecx_setupdatagram(port, &(port->txbuf[idx]), EC_CMD_FPRD, idx,
      *(configlst + slcnt), ECT_REG_ALSTAT, sizeof(ec_alstatust), slstatlst + slcnt);
   sldatapos[slcnt] = EC_HEADERSIZE;
   while(++slcnt < (n - 1))
   {
      sldatapos[slcnt] = ecx_adddatagram(port, &(port->txbuf[idx]), EC_CMD_FPRD, idx, TRUE,
                            *(configlst + slcnt), ECT_REG_ALSTAT, sizeof(ec_alstatust), slstatlst + slcnt);
   }
   if(slcnt < n)
   {
      sldatapos[slcnt] = ecx_adddatagram(port, &(port->txbuf[idx]), EC_CMD_FPRD, idx, FALSE,
                            *(configlst + slcnt), ECT_REG_ALSTAT, sizeof(ec_alstatust), slstatlst + slcnt);
   }
   wkc = ecx_srconfirm(port, idx, timeout);
   if (wkc >= 0)
   {
      for(slcnt = 0 ; slcnt < n ; slcnt++)
      {
         memcpy(slstatlst + slcnt, &(port->rxbuf[idx][sldatapos[slcnt]]), sizeof(ec_alstatust));
      }
   }
   ecx_setbufstat(port, idx, EC_BUF_EMPTY);
   return wkc;
}

/** Read all slave states in ec_slave.
 * @param[in] context = context struct
 * @return lowest state found
 */
int ecx_readstate(ecx_contextt *context)
{
   uint16 slave, fslave, lslave, configadr, lowest, rval;
   ec_alstatust sl[MAX_FPRD_MULTI];
   uint16 slca[MAX_FPRD_MULTI];

   lowest = 0xff;
   context->slavelist[0].ALstatuscode = 0;
   fslave = 1;
   do
   {
      lslave = *(context->slavecount);
      if ((lslave - fslave) >= MAX_FPRD_MULTI)
      {
         lslave = fslave + MAX_FPRD_MULTI - 1;
      }
      for (slave = fslave; slave <= lslave; slave++)
      {
         const ec_alstatust zero = {0, 0, 0};

         configadr = context->slavelist[slave].configadr;
         slca[slave - fslave] = configadr;
         sl[slave - fslave] = zero;
      }
      ecx_FPRD_multi(context, (lslave - fslave) + 1, &(slca[0]), &(sl[0]), EC_TIMEOUTRET3);
      for (slave = fslave; slave <= lslave; slave++)
      {
         configadr = context->slavelist[slave].configadr;
         rval = etohs(sl[slave - fslave].alstatus);
         context->slavelist[slave].ALstatuscode = etohs(sl[slave - fslave].alstatuscode);
         if ((rval & 0xf) < lowest)
         {
            lowest = (rval & 0xf);
         }
         context->slavelist[slave].state = rval;
         context->slavelist[0].ALstatuscode |= context->slavelist[slave].ALstatuscode;
      }
      fslave = lslave + 1;
   } while(lslave < *(context->slavecount));
   context->slavelist[0].state = lowest;

   return lowest;
}

/** Write slave state, if slave = 0 then write to all slaves.
 * The function does not check if the actual state is changed.
 * @param[in]  context        = context struct
 * @param[in] slave    = Slave number, 0 = master
 * @return Workcounter or EC_NOFRAME
 */
int ecx_writestate(ecx_contextt *context, uint16 slave)
{
   int ret;
   uint16 configadr, slstate;

   if (slave == 0)
   {
      slstate = htoes(context->slavelist[slave].state);
      ret = ecx_BWR(context->port, 0, ECT_REG_ALCTL, sizeof(slstate),
	            &slstate, EC_TIMEOUTRET3);
   }
   else
   {
      configadr = context->slavelist[slave].configadr;

      ret = ecx_FPWRw(context->port, configadr, ECT_REG_ALCTL,
	        htoes(context->slavelist[slave].state), EC_TIMEOUTRET3);
   }
   return ret;
}

/** Check actual slave state.
 * This is a blocking function.
 * @param[in] context     = context struct
 * @param[in] slave       = Slave number, 0 = all slaves
 * @param[in] reqstate    = Requested state
 * @param[in] timeout     = Timout value in us
 * @return Requested state, or found state after timeout.
 */
uint16 ecx_statecheck(ecx_contextt *context, uint16 slave, uint16 reqstate, int timeout)
{
   uint16 configadr, state, rval;
   ec_alstatust slstat;
   osal_timert timer;

   if ( slave > *(context->slavecount) )
   {
      return 0;
   }
   osal_timer_start(&timer, timeout);
   configadr = context->slavelist[slave].configadr;
   do
   {
      if (slave < 1)
      {
         rval = 0;
         ecx_BRD(context->port, 0, ECT_REG_ALSTAT, sizeof(rval), &rval , EC_TIMEOUTRET);
         rval = etohs(rval);
      }
      else
      {
         slstat.alstatus = 0;
         slstat.alstatuscode = 0;
         ecx_FPRD(context->port, configadr, ECT_REG_ALSTAT, sizeof(slstat), &slstat, EC_TIMEOUTRET);
         rval = etohs(slstat.alstatus);
         context->slavelist[slave].ALstatuscode = etohs(slstat.alstatuscode);
      }
      state = rval & 0x000f; /* read slave status */
      if (state != reqstate)
      {
         osal_usleep(1000);
      }
   }
   while ((state != reqstate) && (osal_timer_is_expired(&timer) == FALSE));
   context->slavelist[slave].state = rval;

   return state;
}

/** Get index of next mailbox counter value.
 * Used for Mailbox Link Layer.
 * @param[in] cnt     = Mailbox counter value [0..7]
 * @return next mailbox counter value
 */
uint8 ec_nextmbxcnt(uint8 cnt)
{
   cnt++;
   if (cnt > 7)
   {
      cnt = 1; /* wrap around to 1, not 0 */
   }

   return cnt;
}

/** Clear mailbox buffer.
 * @param[out] Mbx     = Mailbox buffer to clear
 */
void ec_clearmbx(ec_mbxbuft *Mbx)
{
    memset(Mbx, 0x00, EC_MAXMBX);
}

/** Check if IN mailbox of slave is empty.
 * @param[in] context  = context struct
 * @param[in] slave    = Slave number
 * @param[in] timeout  = Timeout in us
 * @return >0 is success
 */
int ecx_mbxempty(ecx_contextt *context, uint16 slave, int timeout)
{
   uint16 configadr;
   uint8 SMstat;
   int wkc;
   osal_timert timer;

   osal_timer_start(&timer, timeout);
   configadr = context->slavelist[slave].configadr;
   do
   {
      SMstat = 0;
      wkc = ecx_FPRD(context->port, configadr, ECT_REG_SM0STAT, sizeof(SMstat), &SMstat, EC_TIMEOUTRET);
      SMstat = etohs(SMstat);
      if (((SMstat & 0x08) != 0) && (timeout > EC_LOCALDELAY))
      {
         osal_usleep(EC_LOCALDELAY);
      }
   }
   while (((wkc <= 0) || ((SMstat & 0x08) != 0)) && (osal_timer_is_expired(&timer) == FALSE));

   if ((wkc > 0) && ((SMstat & 0x08) == 0))
   {
      return 1;
   }

   return 0;
}

/** Write IN mailbox to slave.
 * @param[in]  context    = context struct
 * @param[in]  slave      = Slave number
 * @param[out] mbx        = Mailbox data
 * @param[in]  timeout    = Timeout in us
 * @return Work counter (>0 is success)
 */
int ecx_mbxsend(ecx_contextt *context, uint16 slave,ec_mbxbuft *mbx, int timeout)
{
   uint16 mbxwo,mbxl,configadr;
   int wkc;

   wkc = 0;
   configadr = context->slavelist[slave].configadr;
   mbxl = context->slavelist[slave].mbx_l;
   if ((mbxl > 0) && (mbxl <= EC_MAXMBX))
   {
      if (ecx_mbxempty(context, slave, timeout))
      {
         mbxwo = context->slavelist[slave].mbx_wo;
         /* write slave in mailbox */
         wkc = ecx_FPWR(context->port, configadr, mbxwo, mbxl, mbx, EC_TIMEOUTRET3);
      }
      else
      {
         wkc = 0;
      }
   }

   return wkc;
}

/** Read OUT mailbox from slave.
 * Supports Mailbox Link Layer with repeat requests.
 * @param[in]  context    = context struct
 * @param[in]  slave      = Slave number
 * @param[out] mbx        = Mailbox data
 * @param[in]  timeout    = Timeout in us
 * @return Work counter (>0 is success)
 */
int ecx_mbxreceive(ecx_contextt *context, uint16 slave, ec_mbxbuft *mbx, int timeout)
{
   uint16 mbxro,mbxl,configadr;
   int wkc=0;
   int wkc2;
   uint16 SMstat;
   uint8 SMcontr;
   ec_mbxheadert *mbxh;
   ec_emcyt *EMp;
   ec_mbxerrort *MBXEp;

   configadr = context->slavelist[slave].configadr;
   mbxl = context->slavelist[slave].mbx_rl;
   if ((mbxl > 0) && (mbxl <= EC_MAXMBX))
   {
      osal_timert timer;

      osal_timer_start(&timer, timeout);
      wkc = 0;
      do /* wait for read mailbox available */
      {
         SMstat = 0;
         wkc = ecx_FPRD(context->port, configadr, ECT_REG_SM1STAT, sizeof(SMstat), &SMstat, EC_TIMEOUTRET);
         SMstat = etohs(SMstat);
         if (((SMstat & 0x08) == 0) && (timeout > EC_LOCALDELAY))
         {
            osal_usleep(EC_LOCALDELAY);
         }
      }
      while (((wkc <= 0) || ((SMstat & 0x08) == 0)) && (osal_timer_is_expired(&timer) == FALSE));

      if ((wkc > 0) && ((SMstat & 0x08) > 0)) /* read mailbox available ? */
      {
         mbxro = context->slavelist[slave].mbx_ro;
         mbxh = (ec_mbxheadert *)mbx;
         do
         {
            wkc = ecx_FPRD(context->port, configadr, mbxro, mbxl, mbx, EC_TIMEOUTRET); /* get mailbox */
            if ((wkc > 0) && ((mbxh->mbxtype & 0x0f) == 0x00)) /* Mailbox error response? */
            {
               MBXEp = (ec_mbxerrort *)mbx;
               ecx_mbxerror(context, slave, etohs(MBXEp->Detail));
               wkc = 0; /* prevent emergency to cascade up, it is already handled. */
            }
            else if ((wkc > 0) && ((mbxh->mbxtype & 0x0f) == 0x03)) /* CoE response? */
            {
               EMp = (ec_emcyt *)mbx;
               if ((etohs(EMp->CANOpen) >> 12) == 0x01) /* Emergency request? */
               {
                  ecx_mbxemergencyerror(context, slave, etohs(EMp->ErrorCode), EMp->ErrorReg,
                          EMp->bData, etohs(EMp->w1), etohs(EMp->w2));
                  wkc = 0; /* prevent emergency to cascade up, it is already handled. */
               }
            }
            else
            {
               if (wkc <= 0) /* read mailbox lost */
               {
                  SMstat ^= 0x0200; /* toggle repeat request */
                  SMstat = htoes(SMstat);
                  wkc2 = ecx_FPWR(context->port, configadr, ECT_REG_SM1STAT, sizeof(SMstat), &SMstat, EC_TIMEOUTRET);
                  SMstat = etohs(SMstat);
                  do /* wait for toggle ack */
                  {
                     wkc2 = ecx_FPRD(context->port, configadr, ECT_REG_SM1CONTR, sizeof(SMcontr), &SMcontr, EC_TIMEOUTRET);
                   } while (((wkc2 <= 0) || ((SMcontr & 0x02) != (HI_BYTE(SMstat) & 0x02))) && (osal_timer_is_expired(&timer) == FALSE));
                  do /* wait for read mailbox available */
                  {
                     wkc2 = ecx_FPRD(context->port, configadr, ECT_REG_SM1STAT, sizeof(SMstat), &SMstat, EC_TIMEOUTRET);
                     SMstat = etohs(SMstat);
                     if (((SMstat & 0x08) == 0) && (timeout > EC_LOCALDELAY))
                     {
                        osal_usleep(EC_LOCALDELAY);
                     }
                  } while (((wkc2 <= 0) || ((SMstat & 0x08) == 0)) && (osal_timer_is_expired(&timer) == FALSE));
               }
            }
         } while ((wkc <= 0) && (osal_timer_is_expired(&timer) == FALSE)); /* if WKC<=0 repeat */
      }
      else /* no read mailbox available */
      {
          wkc = 0;
      }
   }

   return wkc;
}

/** Dump complete EEPROM data from slave in buffer.
 * @param[in]  context  = context struct
 * @param[in]  slave    = Slave number
 * @param[out] esibuf   = EEPROM data buffer, make sure it is big enough.
 */
void ecx_esidump(ecx_contextt *context, uint16 slave, uint8 *esibuf)
{
   int address, incr;
   uint16 configadr;
   uint64 *p64;
   uint16 *p16;
   uint64 edat;
   uint8 eectl = context->slavelist[slave].eep_pdi;

   ecx_eeprom2master(context, slave); /* set eeprom control to master */
   configadr = context->slavelist[slave].configadr;
   address = ECT_SII_START;
   p16=(uint16*)esibuf;
   if (context->slavelist[slave].eep_8byte)
   {
      incr = 4;
   }
   else
   {
      incr = 2;
   }
   do
   {
      edat = ecx_readeepromFP(context, configadr, address, EC_TIMEOUTEEP);
      p64 = (uint64*)p16;
      *p64 = edat;
      p16 += incr;
      address += incr;
   } while ((address <= (EC_MAXEEPBUF >> 1)) && ((uint32)edat != 0xffffffff));

   if (eectl)
   {
      ecx_eeprom2pdi(context, slave); /* if eeprom control was previously pdi then restore */
   }
}

/** Read EEPROM from slave bypassing cache.
 * @param[in] context   = context struct
 * @param[in] slave     = Slave number
 * @param[in] eeproma   = (WORD) Address in the EEPROM
 * @param[in] timeout   = Timeout in us.
 * @return EEPROM data 32bit
 */
uint32 ecx_readeeprom(ecx_contextt *context, uint16 slave, uint16 eeproma, int timeout)
{
   uint16 configadr;

   ecx_eeprom2master(context, slave); /* set eeprom control to master */
   configadr = context->slavelist[slave].configadr;

   return ((uint32)ecx_readeepromFP(context, configadr, eeproma, timeout));
}

/** Write EEPROM to slave bypassing cache.
 * @param[in] context   = context struct
 * @param[in] slave     = Slave number
 * @param[in] eeproma   = (WORD) Address in the EEPROM
 * @param[in] data      = 16bit data
 * @param[in] timeout   = Timeout in us.
 * @return >0 if OK
 */
int ecx_writeeeprom(ecx_contextt *context, uint16 slave, uint16 eeproma, uint16 data, int timeout)
{
   uint16 configadr;

   ecx_eeprom2master(context, slave); /* set eeprom control to master */
   configadr = context->slavelist[slave].configadr;
   return (ecx_writeeepromFP(context, configadr, eeproma, data, timeout));
}

/** Set eeprom control to master. Only if set to PDI.
 * @param[in] context   = context struct
 * @param[in] slave     = Slave number
 * @return >0 if OK
 */
int ecx_eeprom2master(ecx_contextt *context, uint16 slave)
{
   int wkc = 1, cnt = 0;
   uint16 configadr;
   uint8 eepctl;

   if ( context->slavelist[slave].eep_pdi )
   {
      configadr = context->slavelist[slave].configadr;
      eepctl = 2;
      do
      {
         wkc = ecx_FPWR(context->port, configadr, ECT_REG_EEPCFG, sizeof(eepctl), &eepctl , EC_TIMEOUTRET); /* force Eeprom from PDI */
      }
      while ((wkc <= 0) && (cnt++ < EC_DEFAULTRETRIES));
      eepctl = 0;
      cnt = 0;
      do
      {
         wkc = ecx_FPWR(context->port, configadr, ECT_REG_EEPCFG, sizeof(eepctl), &eepctl , EC_TIMEOUTRET); /* set Eeprom to master */
      }
      while ((wkc <= 0) && (cnt++ < EC_DEFAULTRETRIES));
      context->slavelist[slave].eep_pdi = 0;
   }

   return wkc;
}

/** Set eeprom control to PDI. Only if set to master.
 * @param[in]  context        = context struct
 * @param[in] slave     = Slave number
 * @return >0 if OK
 */
int ecx_eeprom2pdi(ecx_contextt *context, uint16 slave)
{
   int wkc = 1, cnt = 0;
   uint16 configadr;
   uint8 eepctl;

   if ( !context->slavelist[slave].eep_pdi )
   {
      configadr = context->slavelist[slave].configadr;
      eepctl = 1;
      do
      {
         wkc = ecx_FPWR(context->port, configadr, ECT_REG_EEPCFG, sizeof(eepctl), &eepctl , EC_TIMEOUTRET); /* set Eeprom to PDI */
      }
      while ((wkc <= 0) && (cnt++ < EC_DEFAULTRETRIES));
      context->slavelist[slave].eep_pdi = 1;
   }

   return wkc;
}

uint16 ecx_eeprom_waitnotbusyAP(ecx_contextt *context, uint16 aiadr,uint16 *estat, int timeout)
{
   int wkc, cnt = 0, retval = 0;
   osal_timert timer;

   osal_timer_start(&timer, timeout);
   do
   {
      if (cnt++)
      {
         osal_usleep(EC_LOCALDELAY);
      }
      *estat = 0;
      wkc=ecx_APRD(context->port, aiadr, ECT_REG_EEPSTAT, sizeof(*estat), estat, EC_TIMEOUTRET);
      *estat = etohs(*estat);
   }
   while (((wkc <= 0) || ((*estat & EC_ESTAT_BUSY) > 0)) && (osal_timer_is_expired(&timer) == FALSE)); /* wait for eeprom ready */
   if ((*estat & EC_ESTAT_BUSY) == 0)
   {
      retval = 1;
   }

   return retval;
}

/** Read EEPROM from slave bypassing cache. APRD method.
 * @param[in] context     = context struct
 * @param[in] aiadr       = auto increment address of slave
 * @param[in] eeproma     = (WORD) Address in the EEPROM
 * @param[in] timeout     = Timeout in us.
 * @return EEPROM data 64bit or 32bit
 */
uint64 ecx_readeepromAP(ecx_contextt *context, uint16 aiadr, uint16 eeproma, int timeout)
{
   uint16 estat;
   uint32 edat32;
   uint64 edat64;
   ec_eepromt ed;
   int wkc, cnt, nackcnt = 0;

   edat64 = 0;
   edat32 = 0;
   if (ecx_eeprom_waitnotbusyAP(context, aiadr, &estat, timeout))
   {
      if (estat & EC_ESTAT_EMASK) /* error bits are set */
      {
         estat = htoes(EC_ECMD_NOP); /* clear error bits */
         wkc = ecx_APWR(context->port, aiadr, ECT_REG_EEPCTL, sizeof(estat), &estat, EC_TIMEOUTRET3);
      }

      do
      {
         ed.comm = htoes(EC_ECMD_READ);
         ed.addr = htoes(eeproma);
         ed.d2   = 0x0000;
         cnt = 0;
         do
         {
            wkc = ecx_APWR(context->port, aiadr, ECT_REG_EEPCTL, sizeof(ed), &ed, EC_TIMEOUTRET);
         }
         while ((wkc <= 0) && (cnt++ < EC_DEFAULTRETRIES));
         if (wkc)
         {
            osal_usleep(EC_LOCALDELAY);
            estat = 0x0000;
            if (ecx_eeprom_waitnotbusyAP(context, aiadr, &estat, timeout))
            {
               if (estat & EC_ESTAT_NACK)
               {
                  nackcnt++;
                  osal_usleep(EC_LOCALDELAY * 5);
               }
               else
               {
                  nackcnt = 0;
                  if (estat & EC_ESTAT_R64)
                  {
                     cnt = 0;
                     do
                     {
                        wkc = ecx_APRD(context->port, aiadr, ECT_REG_EEPDAT, sizeof(edat64), &edat64, EC_TIMEOUTRET);
                     }
                     while ((wkc <= 0) && (cnt++ < EC_DEFAULTRETRIES));
                  }
                  else
                  {
                     cnt = 0;
                     do
                     {
                        wkc = ecx_APRD(context->port, aiadr, ECT_REG_EEPDAT, sizeof(edat32), &edat32, EC_TIMEOUTRET);
                     }
                     while ((wkc <= 0) && (cnt++ < EC_DEFAULTRETRIES));
                     edat64=(uint64)edat32;
                  }
               }
            }
         }
      }
      while ((nackcnt > 0) && (nackcnt < 3));
   }

   return edat64;
}

/** Write EEPROM to slave bypassing cache. APWR method.
 * @param[in] context   = context struct
 * @param[in] aiadr     = configured address of slave
 * @param[in] eeproma   = (WORD) Address in the EEPROM
 * @param[in] data      = 16bit data
 * @param[in] timeout   = Timeout in us.
 * @return >0 if OK
 */
int ecx_writeeepromAP(ecx_contextt *context, uint16 aiadr, uint16 eeproma, uint16 data, int timeout)
{
   uint16 estat;
   ec_eepromt ed;
   int wkc, rval = 0, cnt = 0, nackcnt = 0;

   if (ecx_eeprom_waitnotbusyAP(context, aiadr, &estat, timeout))
   {
      if (estat & EC_ESTAT_EMASK) /* error bits are set */
      {
         estat = htoes(EC_ECMD_NOP); /* clear error bits */
         wkc = ecx_APWR(context->port, aiadr, ECT_REG_EEPCTL, sizeof(estat), &estat, EC_TIMEOUTRET3);
      }
      do
      {
         cnt = 0;
         do
         {
            wkc = ecx_APWR(context->port, aiadr, ECT_REG_EEPDAT, sizeof(data), &data, EC_TIMEOUTRET);
         }
         while ((wkc <= 0) && (cnt++ < EC_DEFAULTRETRIES));

         ed.comm = EC_ECMD_WRITE;
         ed.addr = eeproma;
         ed.d2   = 0x0000;
         cnt = 0;
         do
         {
            wkc = ecx_APWR(context->port, aiadr, ECT_REG_EEPCTL, sizeof(ed), &ed, EC_TIMEOUTRET);
         }
         while ((wkc <= 0) && (cnt++ < EC_DEFAULTRETRIES));
         if (wkc)
         {
            osal_usleep(EC_LOCALDELAY * 2);
            estat = 0x0000;
            if (ecx_eeprom_waitnotbusyAP(context, aiadr, &estat, timeout))
            {
               if (estat & EC_ESTAT_NACK)
               {
                  nackcnt++;
                  osal_usleep(EC_LOCALDELAY * 5);
               }
               else
               {
                  nackcnt = 0;
                  rval = 1;
               }
            }
         }

      }
      while ((nackcnt > 0) && (nackcnt < 3));
   }

   return rval;
}

uint16 ecx_eeprom_waitnotbusyFP(ecx_contextt *context, uint16 configadr,uint16 *estat, int timeout)
{
   int wkc, cnt = 0, retval = 0;
   osal_timert timer;

   osal_timer_start(&timer, timeout);
   do
   {
      if (cnt++)
      {
         osal_usleep(EC_LOCALDELAY);
      }
      *estat = 0;
      wkc=ecx_FPRD(context->port, configadr, ECT_REG_EEPSTAT, sizeof(*estat), estat, EC_TIMEOUTRET);
      *estat = etohs(*estat);
   }
   while (((wkc <= 0) || ((*estat & EC_ESTAT_BUSY) > 0)) && (osal_timer_is_expired(&timer) == FALSE)); /* wait for eeprom ready */
   if ((*estat & EC_ESTAT_BUSY) == 0)
   {
      retval = 1;
   }

   return retval;
}

/** Read EEPROM from slave bypassing cache. FPRD method.
 * @param[in] context     = context struct
 * @param[in] configadr   = configured address of slave
 * @param[in] eeproma     = (WORD) Address in the EEPROM
 * @param[in] timeout     = Timeout in us.
 * @return EEPROM data 64bit or 32bit
 */
uint64 ecx_readeepromFP(ecx_contextt *context, uint16 configadr, uint16 eeproma, int timeout)
{
   uint16 estat;
   uint32 edat32;
   uint64 edat64;
   ec_eepromt ed;
   int wkc, cnt, nackcnt = 0;

   edat64 = 0;
   edat32 = 0;
   if (ecx_eeprom_waitnotbusyFP(context, configadr, &estat, timeout))
   {
      if (estat & EC_ESTAT_EMASK) /* error bits are set */
      {
         estat = htoes(EC_ECMD_NOP); /* clear error bits */
         wkc=ecx_FPWR(context->port, configadr, ECT_REG_EEPCTL, sizeof(estat), &estat, EC_TIMEOUTRET3);
      }

      do
      {
         ed.comm = htoes(EC_ECMD_READ);
         ed.addr = htoes(eeproma);
         ed.d2   = 0x0000;
         cnt = 0;
         do
         {
            wkc=ecx_FPWR(context->port, configadr, ECT_REG_EEPCTL, sizeof(ed), &ed, EC_TIMEOUTRET);
         }
         while ((wkc <= 0) && (cnt++ < EC_DEFAULTRETRIES));
         if (wkc)
         {
            osal_usleep(EC_LOCALDELAY);
            estat = 0x0000;
            if (ecx_eeprom_waitnotbusyFP(context, configadr, &estat, timeout))
            {
               if (estat & EC_ESTAT_NACK)
               {
                  nackcnt++;
                  osal_usleep(EC_LOCALDELAY * 5);
               }
               else
               {
                  nackcnt = 0;
                  if (estat & EC_ESTAT_R64)
                  {
                     cnt = 0;
                     do
                     {
                        wkc=ecx_FPRD(context->port, configadr, ECT_REG_EEPDAT, sizeof(edat64), &edat64, EC_TIMEOUTRET);
                     }
                     while ((wkc <= 0) && (cnt++ < EC_DEFAULTRETRIES));
                  }
                  else
                  {
                     cnt = 0;
                     do
                     {
                        wkc=ecx_FPRD(context->port, configadr, ECT_REG_EEPDAT, sizeof(edat32), &edat32, EC_TIMEOUTRET);
                     }
                     while ((wkc <= 0) && (cnt++ < EC_DEFAULTRETRIES));
                     edat64=(uint64)edat32;
                  }
               }
            }
         }
      }
      while ((nackcnt > 0) && (nackcnt < 3));
   }

   return edat64;
}

/** Write EEPROM to slave bypassing cache. FPWR method.
 * @param[in]  context        = context struct
 * @param[in] configadr   = configured address of slave
 * @param[in] eeproma     = (WORD) Address in the EEPROM
 * @param[in] data        = 16bit data
 * @param[in] timeout     = Timeout in us.
 * @return >0 if OK
 */
int ecx_writeeepromFP(ecx_contextt *context, uint16 configadr, uint16 eeproma, uint16 data, int timeout)
{
   uint16 estat;
   ec_eepromt ed;
   int wkc, rval = 0, cnt = 0, nackcnt = 0;

   if (ecx_eeprom_waitnotbusyFP(context, configadr, &estat, timeout))
   {
      if (estat & EC_ESTAT_EMASK) /* error bits are set */
      {
         estat = htoes(EC_ECMD_NOP); /* clear error bits */
         wkc = ecx_FPWR(context->port, configadr, ECT_REG_EEPCTL, sizeof(estat), &estat, EC_TIMEOUTRET3);
      }
      do
      {
         cnt = 0;
         do
         {
            wkc = ecx_FPWR(context->port, configadr, ECT_REG_EEPDAT, sizeof(data), &data, EC_TIMEOUTRET);
         }
         while ((wkc <= 0) && (cnt++ < EC_DEFAULTRETRIES));
         ed.comm = EC_ECMD_WRITE;
         ed.addr = eeproma;
         ed.d2   = 0x0000;
         cnt = 0;
         do
         {
            wkc = ecx_FPWR(context->port, configadr, ECT_REG_EEPCTL, sizeof(ed), &ed, EC_TIMEOUTRET);
         }
         while ((wkc <= 0) && (cnt++ < EC_DEFAULTRETRIES));
         if (wkc)
         {
            osal_usleep(EC_LOCALDELAY * 2);
            estat = 0x0000;
            if (ecx_eeprom_waitnotbusyFP(context, configadr, &estat, timeout))
            {
               if (estat & EC_ESTAT_NACK)
               {
                  nackcnt++;
                  osal_usleep(EC_LOCALDELAY * 5);
               }
               else
               {
                  nackcnt = 0;
                  rval = 1;
               }
            }
         }
      }
      while ((nackcnt > 0) && (nackcnt < 3));
   }

   return rval;
}

/** Read EEPROM from slave bypassing cache.
 * Parallel read step 1, make request to slave.
 * @param[in] context     = context struct
 * @param[in] slave       = Slave number
 * @param[in] eeproma     = (WORD) Address in the EEPROM
 */
void ecx_readeeprom1(ecx_contextt *context, uint16 slave, uint16 eeproma)
{
   uint16 configadr, estat;
   ec_eepromt ed;
   int wkc, cnt = 0;

   ecx_eeprom2master(context, slave); /* set eeprom control to master */
   configadr = context->slavelist[slave].configadr;
   if (ecx_eeprom_waitnotbusyFP(context, configadr, &estat, EC_TIMEOUTEEP))
   {
      if (estat & EC_ESTAT_EMASK) /* error bits are set */
      {
         estat = htoes(EC_ECMD_NOP); /* clear error bits */
         wkc = ecx_FPWR(context->port, configadr, ECT_REG_EEPCTL, sizeof(estat), &estat, EC_TIMEOUTRET3);
      }
      ed.comm = htoes(EC_ECMD_READ);
      ed.addr = htoes(eeproma);
      ed.d2   = 0x0000;
      do
      {
         wkc = ecx_FPWR(context->port, configadr, ECT_REG_EEPCTL, sizeof(ed), &ed, EC_TIMEOUTRET);
      }
      while ((wkc <= 0) && (cnt++ < EC_DEFAULTRETRIES));
   }
}

/** Read EEPROM from slave bypassing cache.
 * Parallel read step 2, actual read from slave.
 * @param[in]  context        = context struct
 * @param[in] slave       = Slave number
 * @param[in] timeout     = Timeout in us.
 * @return EEPROM data 32bit
 */
uint32 ecx_readeeprom2(ecx_contextt *context, uint16 slave, int timeout)
{
   uint16 estat, configadr;
   uint32 edat;
   int wkc, cnt = 0;

   configadr = context->slavelist[slave].configadr;
   edat = 0;
   estat = 0x0000;
   if (ecx_eeprom_waitnotbusyFP(context, configadr, &estat, timeout))
   {
      do
      {
          wkc = ecx_FPRD(context->port, configadr, ECT_REG_EEPDAT, sizeof(edat), &edat, EC_TIMEOUTRET);
      }
      while ((wkc <= 0) && (cnt++ < EC_DEFAULTRETRIES));
   }

   return edat;
}

/** Push index of segmented LRD/LWR/LRW combination.
 * @param[in]  context        = context struct
 * @param[in] idx         = Used datagram index.
 * @param[in] data        = Pointer to process data segment.
 * @param[in] length      = Length of data segment in bytes.
 */
static void ecx_pushindex(ecx_contextt *context, uint8 idx, void *data, uint16 length)
{
   if(context->idxstack->pushed < EC_MAXBUF)
   {
      context->idxstack->idx[context->idxstack->pushed] = idx;
      context->idxstack->data[context->idxstack->pushed] = data;
      context->idxstack->length[context->idxstack->pushed] = length;
      context->idxstack->pushed++;
   }
}

/** Pull index of segmented LRD/LWR/LRW combination.
 * @param[in]  context        = context struct
 * @return Stack location, -1 if stack is empty.
 */
static int ecx_pullindex(ecx_contextt *context)
{
   int rval = -1;
   if(context->idxstack->pulled < context->idxstack->pushed)
   {
      rval = context->idxstack->pulled;
      context->idxstack->pulled++;
   }

   return rval;
}

/** 
 * Clear the idx stack.
 * 
 * @param context           = context struct
 */
static void ecx_clearindex(ecx_contextt *context)  {

   context->idxstack->pushed = 0;
   context->idxstack->pulled = 0;

}

/** Transmit processdata to slaves.
 * Uses LRW, or LRD/LWR if LRW is not allowed (blockLRW).
 * Both the input and output processdata are transmitted.
 * The outputs with the actual data, the inputs have a placeholder.
 * The inputs are gathered with the receive processdata function.
 * In contrast to the base LRW function this function is non-blocking.
 * If the processdata does not fit in one datagram, multiple are used.
 * In order to recombine the slave response, a stack is used.
 * @param[in]  context        = context struct
 * @param[in]  group          = group number
 * @return >0 if processdata is transmitted.
 */
int ecx_send_processdata_group(ecx_contextt *context, uint8 group)
{
   uint32 LogAdr;
   uint16 w1, w2;
   int length, sublength;
   uint8 idx;
   int wkc;
   uint8* data;
   boolean first=FALSE;
   uint16 currentsegment = 0;

   wkc = 0;
   if(context->grouplist[group].hasdc)
   {
      first = TRUE;
   }
   length = context->grouplist[group].Obytes + context->grouplist[group].Ibytes;
   LogAdr = context->grouplist[group].logstartaddr;
   if (length)
   {

      wkc = 1;
      /* LRW blocked by one or more slaves ? */
      if (context->grouplist[group].blockLRW)
      {
         /* if inputs available generate LRD */
         if(context->grouplist[group].Ibytes)
         {
            currentsegment = context->grouplist[group].Isegment;
            data = context->grouplist[group].inputs;
            length = context->grouplist[group].Ibytes;
            LogAdr += context->grouplist[group].Obytes;
            /* segment transfer if needed */
            do
            {
               if(currentsegment == context->grouplist[group].Isegment)
               {
                  sublength = context->grouplist[group].IOsegment[currentsegment++] - context->grouplist[group].Ioffset;
               }
               else
               {
                  sublength = context->grouplist[group].IOsegment[currentsegment++];
               }
               /* get new index */
               idx = ecx_getindex(context->port);
               w1 = LO_WORD(LogAdr);
               w2 = HI_WORD(LogAdr);
               ecx_setupdatagram(context->port, &(context->port->txbuf[idx]), EC_CMD_LRD, idx, w1, w2, sublength, data);
               if(first)
               {
                  context->DCl = sublength;
                  /* FPRMW in second datagram */
                  context->DCtO = ecx_adddatagram(context->port, &(context->port->txbuf[idx]), EC_CMD_FRMW, idx, FALSE,
                                           context->slavelist[context->grouplist[group].DCnext].configadr,
                                           ECT_REG_DCSYSTIME, sizeof(int64), context->DCtime);
                  first = FALSE;
               }
               /* send frame */
               ecx_outframe_red(context->port, idx);
               /* push index and data pointer on stack */
               ecx_pushindex(context, idx, data, sublength);
               length -= sublength;
               LogAdr += sublength;
               data += sublength;
            } while (length && (currentsegment < context->grouplist[group].nsegments));
         }
         /* if outputs available generate LWR */
         if(context->grouplist[group].Obytes)
         {
            data = context->grouplist[group].outputs;
            length = context->grouplist[group].Obytes;
            LogAdr = context->grouplist[group].logstartaddr;
            currentsegment = 0;
            /* segment transfer if needed */
            do
            {
               sublength = context->grouplist[group].IOsegment[currentsegment++];
               if((length - sublength) < 0)
               {
                  sublength = length;
               }
               /* get new index */
               idx = ecx_getindex(context->port);
               w1 = LO_WORD(LogAdr);
               w2 = HI_WORD(LogAdr);
               ecx_setupdatagram(context->port, &(context->port->txbuf[idx]), EC_CMD_LWR, idx, w1, w2, sublength, data);
               if(first)
               {
                  context->DCl = sublength;
                  /* FPRMW in second datagram */
                  context->DCtO = ecx_adddatagram(context->port, &(context->port->txbuf[idx]), EC_CMD_FRMW, idx, FALSE,
                                           context->slavelist[context->grouplist[group].DCnext].configadr,
                                           ECT_REG_DCSYSTIME, sizeof(int64), context->DCtime);
                  first = FALSE;
               }
               /* send frame */
               ecx_outframe_red(context->port, idx);
               /* push index and data pointer on stack */
               ecx_pushindex(context, idx, data, sublength);
               length -= sublength;
               LogAdr += sublength;
               data += sublength;
            } while (length && (currentsegment < context->grouplist[group].nsegments));
         }
      }
      /* LRW can be used */
      else
      {
         if (context->grouplist[group].Obytes)
         {
            data = context->grouplist[group].outputs;
         }
         else
         {
            data = context->grouplist[group].inputs;
         }
         /* segment transfer if needed */
         do
         {
            sublength = context->grouplist[group].IOsegment[currentsegment++];
            /* get new index */
            idx = ecx_getindex(context->port);
            w1 = LO_WORD(LogAdr);
            w2 = HI_WORD(LogAdr);
            ecx_setupdatagram(context->port, &(context->port->txbuf[idx]), EC_CMD_LRW, idx, w1, w2, sublength, data);
            if(first)
            {
               context->DCl = sublength;
               /* FPRMW in second datagram */
               context->DCtO = ecx_adddatagram(context->port, &(context->port->txbuf[idx]), EC_CMD_FRMW, idx, FALSE,
                                        context->slavelist[context->grouplist[group].DCnext].configadr,
                                        ECT_REG_DCSYSTIME, sizeof(int64), context->DCtime);
               first = FALSE;
            }
            /* send frame */
            ecx_outframe_red(context->port, idx);
            /* push index and data pointer on stack */
            ecx_pushindex(context, idx, data, sublength);
            length -= sublength;
            LogAdr += sublength;
            data += sublength;
         } while (length && (currentsegment < context->grouplist[group].nsegments));
      }
   }

   return wkc;
}

/** Receive processdata from slaves.
 * Second part from ec_send_processdata().
 * Received datagrams are recombined with the processdata with help from the stack.
 * If a datagram contains input processdata it copies it to the processdata structure.
 * @param[in]  context        = context struct
 * @param[in]  group          = group number
 * @param[in]  timeout        = Timeout in us.
 * @return Work counter.
 */
int ecx_receive_processdata_group(ecx_contextt *context, uint8 group, int timeout)
{
   int pos, idx;
   int wkc = 0, wkc2;
   uint16 le_wkc = 0;
   int valid_wkc = 0;
   int64 le_DCtime;
   boolean first = FALSE;

   if(context->grouplist[group].hasdc)
   {
      first = TRUE;
   }
   /* get first index */
   pos = ecx_pullindex(context);
   /* read the same number of frames as send */
   while (pos >= 0)
   {
      idx = context->idxstack->idx[pos];
      wkc2 = ecx_waitinframe(context->port, context->idxstack->idx[pos], timeout);
      /* check if there is input data in frame */
      if (wkc2 > EC_NOFRAME)
      {
         if((context->port->rxbuf[idx][EC_CMDOFFSET]==EC_CMD_LRD) || (context->port->rxbuf[idx][EC_CMDOFFSET]==EC_CMD_LRW))
         {
            if(first)
            {
               memcpy(context->idxstack->data[pos], &(context->port->rxbuf[idx][EC_HEADERSIZE]), context->DCl);
               memcpy(&le_wkc, &(context->port->rxbuf[idx][EC_HEADERSIZE + context->DCl]), EC_WKCSIZE);
               wkc = etohs(le_wkc);
               memcpy(&le_DCtime, &(context->port->rxbuf[idx][context->DCtO]), sizeof(le_DCtime));
               *(context->DCtime) = etohll(le_DCtime);
               first = FALSE;
            }
            else
            {
               /* copy input data back to process data buffer */
               memcpy(context->idxstack->data[pos], &(context->port->rxbuf[idx][EC_HEADERSIZE]), context->idxstack->length[pos]);
               wkc += wkc2;
            }
            valid_wkc = 1;
         }
         else if(context->port->rxbuf[idx][EC_CMDOFFSET]==EC_CMD_LWR)
         {
            if(first)
            {
               memcpy(&le_wkc, &(context->port->rxbuf[idx][EC_HEADERSIZE + context->DCl]), EC_WKCSIZE);
               /* output WKC counts 2 times when using LRW, emulate the same for LWR */
               wkc = etohs(le_wkc) * 2;
               memcpy(&le_DCtime, &(context->port->rxbuf[idx][context->DCtO]), sizeof(le_DCtime));
               *(context->DCtime) = etohll(le_DCtime);
               first = FALSE;
            }
            else
            {
               /* output WKC counts 2 times when using LRW, emulate the same for LWR */
               wkc += wkc2 * 2;
            }
            valid_wkc = 1;
         }
      }
      /* release buffer */
      ecx_setbufstat(context->port, idx, EC_BUF_EMPTY);
      /* get next index */
      pos = ecx_pullindex(context);
   }

   ecx_clearindex(context);

   /* if no frames has arrived */
   if (valid_wkc == 0)
   {
      return EC_NOFRAME;
   }
   return wkc;
}


int ecx_send_processdata(ecx_contextt *context)
{
   return ecx_send_processdata_group(context, 0);
}

int ecx_receive_processdata(ecx_contextt *context, int timeout)
{
   return ecx_receive_processdata_group(context, 0, timeout);
}

#ifdef EC_VER1
void ec_pusherror(const ec_errort *Ec)
{
   ecx_pusherror(&ecx_context, Ec);
}

boolean ec_poperror(ec_errort *Ec)
{
   return ecx_poperror(&ecx_context, Ec);
}

boolean ec_iserror(void)
{
   return ecx_iserror(&ecx_context);
}

void ec_packeterror(uint16 Slave, uint16 Index, uint8 SubIdx, uint16 ErrorCode)
{
   ecx_packeterror(&ecx_context, Slave, Index, SubIdx, ErrorCode);
}

/** Initialise lib in single NIC mode
 * @param[in] ifname   = Dev name, f.e. "eth0"
 * @return >0 if OK
 * @see ecx_init
 */
int ec_init(const char * ifname)
{
   return ecx_init(&ecx_context, ifname);
}

/** Initialise lib in redundant NIC mode
 * @param[in]  ifname   = Primary Dev name, f.e. "eth0"
 * @param[in]  if2name  = Secondary Dev name, f.e. "eth1"
 * @return >0 if OK
 * @see ecx_init_redundant
 */
int ec_init_redundant(const char *ifname, char *if2name)
{
   return ecx_init_redundant (&ecx_context, &ecx_redport, ifname, if2name);
}

/** Close lib.
 * @see ecx_close
 */
void ec_close(void)
{
   ecx_close(&ecx_context);
};

/** Read one byte from slave EEPROM via cache.
 *  If the cache location is empty then a read request is made to the slave.
 *  Depending on the slave capabillities the request is 4 or 8 bytes.
 *  @param[in] slave   = slave number
 *  @param[in] address = eeprom address in bytes (slave uses words)
 *  @return requested byte, if not available then 0xff
 * @see ecx_siigetbyte
 */
uint8 ec_siigetbyte(uint16 slave, uint16 address)
{
   return ecx_siigetbyte (&ecx_context, slave, address);
}

/** Find SII section header in slave EEPROM.
 *  @param[in] slave   = slave number
 *  @param[in] cat     = section category
 *  @return byte address of section at section length entry, if not available then 0
 *  @see ecx_siifind
 */
int16 ec_siifind(uint16 slave, uint16 cat)
{
   return ecx_siifind (&ecx_context, slave, cat);
}

/** Get string from SII string section in slave EEPROM.
 *  @param[out] str    = requested string, 0x00 if not found
 *  @param[in]  slave  = slave number
 *  @param[in]  Sn     = string number
 *  @see ecx_siistring
 */
void ec_siistring(char *str, uint16 slave, uint16 Sn)
{
   ecx_siistring(&ecx_context, str, slave, Sn);
}

/** Get FMMU data from SII FMMU section in slave EEPROM.
 *  @param[in]  slave  = slave number
 *  @param[out] FMMU   = FMMU struct from SII, max. 4 FMMU's
 *  @return number of FMMU's defined in section
 *  @see ecx_siiFMMU
 */
uint16 ec_siiFMMU(uint16 slave, ec_eepromFMMUt* FMMU)
{
   return ecx_siiFMMU (&ecx_context, slave, FMMU);
}

/** Get SM data from SII SM section in slave EEPROM.
 *  @param[in]  slave   = slave number
 *  @param[out] SM      = first SM struct from SII
 *  @return number of SM's defined in section
 *  @see ecx_siiSM
 */
uint16 ec_siiSM(uint16 slave, ec_eepromSMt* SM)
{
   return ecx_siiSM (&ecx_context, slave, SM);
}

/** Get next SM data from SII SM section in slave EEPROM.
 *  @param[in]  slave  = slave number
 *  @param[out] SM     = first SM struct from SII
 *  @param[in]  n      = SM number
 *  @return >0 if OK
 *  @see ecx_siiSMnext
 */
uint16 ec_siiSMnext(uint16 slave, ec_eepromSMt* SM, uint16 n)
{
   return ecx_siiSMnext (&ecx_context, slave, SM, n);
}

/** Get PDO data from SII PDO section in slave EEPROM.
 *  @param[in]  slave  = slave number
 *  @param[out] PDO    = PDO struct from SII
 *  @param[in]  t      = 0=RXPDO 1=TXPDO
 *  @return mapping size in bits of PDO
 *  @see ecx_siiPDO
 */
int ec_siiPDO(uint16 slave, ec_eepromPDOt* PDO, uint8 t)
{
   return ecx_siiPDO (&ecx_context, slave, PDO, t);
}

/** Read all slave states in ec_slave.
 * @return lowest state found
 * @see ecx_readstate
 */
int ec_readstate(void)
{
   return ecx_readstate (&ecx_context);
}

/** Write slave state, if slave = 0 then write to all slaves.
 * The function does not check if the actual state is changed.
 * @param[in] slave = Slave number, 0 = master
 * @return 0
 * @see ecx_writestate
 */
int ec_writestate(uint16 slave)
{
   return ecx_writestate(&ecx_context, slave);
}

/** Check actual slave state.
 * This is a blocking function.
 * @param[in] slave       = Slave number, 0 = all slaves
 * @param[in] reqstate    = Requested state
 * @param[in] timeout     = Timout value in us
 * @return Requested state, or found state after timeout.
 * @see ecx_statecheck
 */
uint16 ec_statecheck(uint16 slave, uint16 reqstate, int timeout)
{
   return ecx_statecheck (&ecx_context, slave, reqstate, timeout);
}

/** Check if IN mailbox of slave is empty.
 * @param[in] slave    = Slave number
 * @param[in] timeout  = Timeout in us
 * @return >0 is success
 * @see ecx_mbxempty
 */
int ec_mbxempty(uint16 slave, int timeout)
{
   return ecx_mbxempty (&ecx_context, slave, timeout);
}

/** Write IN mailbox to slave.
 * @param[in]  slave      = Slave number
 * @param[out] mbx        = Mailbox data
 * @param[in]  timeout    = Timeout in us
 * @return Work counter (>0 is success)
 * @see ecx_mbxsend
 */
int ec_mbxsend(uint16 slave,ec_mbxbuft *mbx, int timeout)
{
   return ecx_mbxsend (&ecx_context, slave, mbx, timeout);
}

/** Read OUT mailbox from slave.
 * Supports Mailbox Link Layer with repeat requests.
 * @param[in]  slave      = Slave number
 * @param[out] mbx        = Mailbox data
 * @param[in]  timeout    = Timeout in us
 * @return Work counter (>0 is success)
 * @see ecx_mbxreceive
 */
int ec_mbxreceive(uint16 slave, ec_mbxbuft *mbx, int timeout)
{
   return ecx_mbxreceive (&ecx_context, slave, mbx, timeout);
}

/** Dump complete EEPROM data from slave in buffer.
 * @param[in]  slave    = Slave number
 * @param[out] esibuf   = EEPROM data buffer, make sure it is big enough.
 * @see ecx_esidump
 */
void ec_esidump(uint16 slave, uint8 *esibuf)
{
   ecx_esidump (&ecx_context, slave, esibuf);
}

/** Read EEPROM from slave bypassing cache.
 * @param[in] slave     = Slave number
 * @param[in] eeproma   = (WORD) Address in the EEPROM
 * @param[in] timeout   = Timeout in us.
 * @return EEPROM data 32bit
 * @see ecx_readeeprom
 */
uint32 ec_readeeprom(uint16 slave, uint16 eeproma, int timeout)
{
   return ecx_readeeprom (&ecx_context, slave, eeproma, timeout);
}

/** Write EEPROM to slave bypassing cache.
 * @param[in] slave     = Slave number
 * @param[in] eeproma   = (WORD) Address in the EEPROM
 * @param[in] data      = 16bit data
 * @param[in] timeout   = Timeout in us.
 * @return >0 if OK
 * @see ecx_writeeeprom
 */
int ec_writeeeprom(uint16 slave, uint16 eeproma, uint16 data, int timeout)
{
   return ecx_writeeeprom (&ecx_context, slave, eeproma, data, timeout);
}

/** Set eeprom control to master. Only if set to PDI.
 * @param[in] slave = Slave number
 * @return >0 if OK
 * @see ecx_eeprom2master
 */
int ec_eeprom2master(uint16 slave)
{
   return ecx_eeprom2master(&ecx_context, slave);
}

int ec_eeprom2pdi(uint16 slave)
{
   return ecx_eeprom2pdi(&ecx_context, slave);
}

uint16 ec_eeprom_waitnotbusyAP(uint16 aiadr,uint16 *estat, int timeout)
{
   return ecx_eeprom_waitnotbusyAP (&ecx_context, aiadr, estat, timeout);
}

/** Read EEPROM from slave bypassing cache. APRD method.
 * @param[in] aiadr       = auto increment address of slave
 * @param[in] eeproma     = (WORD) Address in the EEPROM
 * @param[in] timeout     = Timeout in us.
 * @return EEPROM data 64bit or 32bit
 */
uint64 ec_readeepromAP(uint16 aiadr, uint16 eeproma, int timeout)
{
   return ecx_readeepromAP (&ecx_context, aiadr, eeproma, timeout);
}

/** Write EEPROM to slave bypassing cache. APWR method.
 * @param[in] aiadr     = configured address of slave
 * @param[in] eeproma   = (WORD) Address in the EEPROM
 * @param[in] data      = 16bit data
 * @param[in] timeout   = Timeout in us.
 * @return >0 if OK
 * @see ecx_writeeepromAP
 */
int ec_writeeepromAP(uint16 aiadr, uint16 eeproma, uint16 data, int timeout)
{
   return ecx_writeeepromAP (&ecx_context, aiadr, eeproma, data, timeout);
}

uint16 ec_eeprom_waitnotbusyFP(uint16 configadr,uint16 *estat, int timeout)
{
   return ecx_eeprom_waitnotbusyFP (&ecx_context, configadr, estat, timeout);
}

/** Read EEPROM from slave bypassing cache. FPRD method.
 * @param[in] configadr   = configured address of slave
 * @param[in] eeproma     = (WORD) Address in the EEPROM
 * @param[in] timeout     = Timeout in us.
 * @return EEPROM data 64bit or 32bit
 * @see ecx_readeepromFP
 */
uint64 ec_readeepromFP(uint16 configadr, uint16 eeproma, int timeout)
{
   return ecx_readeepromFP (&ecx_context, configadr, eeproma, timeout);
}

/** Write EEPROM to slave bypassing cache. FPWR method.
 * @param[in] configadr   = configured address of slave
 * @param[in] eeproma     = (WORD) Address in the EEPROM
 * @param[in] data        = 16bit data
 * @param[in] timeout     = Timeout in us.
 * @return >0 if OK
 * @see ecx_writeeepromFP
 */
int ec_writeeepromFP(uint16 configadr, uint16 eeproma, uint16 data, int timeout)
{
   return ecx_writeeepromFP (&ecx_context, configadr, eeproma, data, timeout);
}

/** Read EEPROM from slave bypassing cache.
 * Parallel read step 1, make request to slave.
 * @param[in] slave       = Slave number
 * @param[in] eeproma     = (WORD) Address in the EEPROM
 * @see ecx_readeeprom1
 */
void ec_readeeprom1(uint16 slave, uint16 eeproma)
{
   ecx_readeeprom1 (&ecx_context, slave, eeproma);
}

/** Read EEPROM from slave bypassing cache.
 * Parallel read step 2, actual read from slave.
 * @param[in] slave       = Slave number
 * @param[in] timeout     = Timeout in us.
 * @return EEPROM data 32bit
 * @see ecx_readeeprom2
 */
uint32 ec_readeeprom2(uint16 slave, int timeout)
{
   return ecx_readeeprom2 (&ecx_context, slave, timeout);
}

/** Transmit processdata to slaves.
 * Uses LRW, or LRD/LWR if LRW is not allowed (blockLRW).
 * Both the input and output processdata are transmitted.
 * The outputs with the actual data, the inputs have a placeholder.
 * The inputs are gathered with the receive processdata function.
 * In contrast to the base LRW function this function is non-blocking.
 * If the processdata does not fit in one datagram, multiple are used.
 * In order to recombine the slave response, a stack is used.
 * @param[in]  group          = group number
 * @return >0 if processdata is transmitted.
 * @see ecx_send_processdata_group
 */
int ec_send_processdata_group(uint8 group)
{
   return ecx_send_processdata_group (&ecx_context, group);
}

/** Receive processdata from slaves.
 * Second part from ec_send_processdata().
 * Received datagrams are recombined with the processdata with help from the stack.
 * If a datagram contains input processdata it copies it to the processdata structure.
 * @param[in]  group          = group number
 * @param[in]  timeout        = Timeout in us.
 * @return Work counter.
 * @see ecx_receive_processdata_group
 */
int ec_receive_processdata_group(uint8 group, int timeout)
{
   return ecx_receive_processdata_group (&ecx_context, group, timeout);
}

int ec_send_processdata(void)
{
   return ec_send_processdata_group(0);
}

int ec_receive_processdata(int timeout)
{
   return ec_receive_processdata_group(0, timeout);
}
#endif
