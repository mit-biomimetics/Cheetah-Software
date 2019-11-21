/*
 * Simple Open EtherCAT Master Library
 *
 * File    : ethercatconfig.c
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
 * Configuration module for EtherCAT master.
 *
 * After successful initialisation with ec_init() or ec_init_redundant()
 * the slaves can be auto configured with this module.
 */

#include <stdio.h>
#include <string.h>
#include "osal.h"
#include "oshw.h"
#include "ethercattype.h"
#include "ethercatbase.h"
#include "ethercatmain.h"
#include "ethercatcoe.h"
#include "ethercatsoe.h"
#include "ethercatconfig.h"

// define if debug printf is needed
//#define EC_DEBUG

#ifdef EC_DEBUG
#define EC_PRINT printf
#else
#define EC_PRINT(...) do {} while (0)
#endif

typedef struct
{
   int thread_n;
   int running;
   ecx_contextt *context;
   uint16 slave;
} ecx_mapt_t;

ecx_mapt_t ecx_mapt[EC_MAX_MAPT];
OSAL_THREAD_HANDLE ecx_threadh[EC_MAX_MAPT];

#ifdef EC_VER1
/** Slave configuration structure */
typedef const struct
{
   /** Manufacturer code of slave */
   uint32           man;
   /** ID of slave */
   uint32           id;
   /** Readable name */
   char             name[EC_MAXNAME + 1];
   /** Data type */
   uint8            Dtype;
   /** Input bits */
   uint16            Ibits;
   /** Output bits */
   uint16           Obits;
   /** SyncManager 2 address */
   uint16           SM2a;
   /** SyncManager 2 flags */
   uint32           SM2f;
   /** SyncManager 3 address */
   uint16           SM3a;
   /** SyncManager 3 flags */
   uint32           SM3f;
   /** FMMU 0 activation */
   uint8            FM0ac;
   /** FMMU 1 activation */
   uint8            FM1ac;
} ec_configlist_t;

#include "ethercatconfiglist.h"
#endif

/** standard SM0 flags configuration for mailbox slaves */
#define EC_DEFAULTMBXSM0  0x00010026
/** standard SM1 flags configuration for mailbox slaves */
#define EC_DEFAULTMBXSM1  0x00010022
/** standard SM0 flags configuration for digital output slaves */
#define EC_DEFAULTDOSM0   0x00010044

#ifdef EC_VER1
/** Find slave in standard configuration list ec_configlist[]
 *
 * @param[in] man      = manufacturer
 * @param[in] id       = ID
 * @return index in ec_configlist[] when found, otherwise 0
 */
int ec_findconfig( uint32 man, uint32 id)
{
   int i = 0;

   do
   {
      i++;
   } while ( (ec_configlist[i].man != EC_CONFIGEND) &&
           ((ec_configlist[i].man != man) || (ec_configlist[i].id != id)) );
   if (ec_configlist[i].man == EC_CONFIGEND)
   {
      i = 0;
   }
   return i;
}
#endif

void ecx_init_context(ecx_contextt *context)
{
   int lp;
   *(context->slavecount) = 0;
   /* clean ec_slave array */
   memset(context->slavelist, 0x00, sizeof(ec_slavet) * context->maxslave);
   memset(context->grouplist, 0x00, sizeof(ec_groupt) * context->maxgroup);
   /* clear slave eeprom cache, does not actually read any eeprom */
   ecx_siigetbyte(context, 0, EC_MAXEEPBUF);
   for(lp = 0; lp < context->maxgroup; lp++)
   {
      context->grouplist[lp].logstartaddr = lp << 16; /* default start address per group entry */
   }
}

int ecx_detect_slaves(ecx_contextt *context)
{
   uint8  b;
   uint16 w;
   int    wkc;

   /* make special pre-init register writes to enable MAC[1] local administered bit *
    * setting for old netX100 slaves */
   b = 0x00;
   ecx_BWR(context->port, 0x0000, ECT_REG_DLALIAS, sizeof(b), &b, EC_TIMEOUTRET3);     /* Ignore Alias register */
   b = EC_STATE_INIT | EC_STATE_ACK;
   ecx_BWR(context->port, 0x0000, ECT_REG_ALCTL, sizeof(b), &b, EC_TIMEOUTRET3);       /* Reset all slaves to Init */
   /* netX100 should now be happy */
   ecx_BWR(context->port, 0x0000, ECT_REG_ALCTL, sizeof(b), &b, EC_TIMEOUTRET3);       /* Reset all slaves to Init */
   wkc = ecx_BRD(context->port, 0x0000, ECT_REG_TYPE, sizeof(w), &w, EC_TIMEOUTSAFE);  /* detect number of slaves */
   if (wkc > 0)
   {
      /* this is strictly "less than" since the master is "slave 0" */
      if (wkc < EC_MAXSLAVE)
      {
         *(context->slavecount) = wkc;
      }
      else
      {
         EC_PRINT("Error: too many slaves on network: num_slaves=%d, EC_MAXSLAVE=%d\n",
               wkc, EC_MAXSLAVE);
         return -2;
      }
   }
   return wkc;
}

static void ecx_set_slaves_to_default(ecx_contextt *context)
{
   uint8 b;
   uint16 w;
   uint8 zbuf[64];
   memset(&zbuf, 0x00, sizeof(zbuf));
   b = 0x00;
   ecx_BWR(context->port, 0x0000, ECT_REG_DLPORT      , sizeof(b) , &b, EC_TIMEOUTRET3);     /* deact loop manual */
   w = htoes(0x0004);
   ecx_BWR(context->port, 0x0000, ECT_REG_IRQMASK     , sizeof(w) , &w, EC_TIMEOUTRET3);     /* set IRQ mask */
   ecx_BWR(context->port, 0x0000, ECT_REG_RXERR       , 8         , &zbuf, EC_TIMEOUTRET3);  /* reset CRC counters */
   ecx_BWR(context->port, 0x0000, ECT_REG_FMMU0       , 16 * 3    , &zbuf, EC_TIMEOUTRET3);  /* reset FMMU's */
   ecx_BWR(context->port, 0x0000, ECT_REG_SM0         , 8 * 4     , &zbuf, EC_TIMEOUTRET3);  /* reset SyncM */
   ecx_BWR(context->port, 0x0000, ECT_REG_DCSYSTIME   , 4         , &zbuf, EC_TIMEOUTRET3);  /* reset system time+ofs */
   w = htoes(0x1000);
   ecx_BWR(context->port, 0x0000, ECT_REG_DCSPEEDCNT  , sizeof(w) , &w, EC_TIMEOUTRET3);     /* DC speedstart */
   w = htoes(0x0c00);
   ecx_BWR(context->port, 0x0000, ECT_REG_DCTIMEFILT  , sizeof(w) , &w, EC_TIMEOUTRET3);     /* DC filt expr */
   b = 0x00;
   ecx_BWR(context->port, 0x0000, ECT_REG_DLALIAS     , sizeof(b) , &b, EC_TIMEOUTRET3);     /* Ignore Alias register */
   b = EC_STATE_INIT | EC_STATE_ACK;
   ecx_BWR(context->port, 0x0000, ECT_REG_ALCTL       , sizeof(b) , &b, EC_TIMEOUTRET3);     /* Reset all slaves to Init */
   b = 2;
   ecx_BWR(context->port, 0x0000, ECT_REG_EEPCFG      , sizeof(b) , &b, EC_TIMEOUTRET3);     /* force Eeprom from PDI */
   b = 0;
   ecx_BWR(context->port, 0x0000, ECT_REG_EEPCFG      , sizeof(b) , &b, EC_TIMEOUTRET3);     /* set Eeprom to master */
}

#ifdef EC_VER1
static int ecx_config_from_table(ecx_contextt *context, uint16 slave)
{
   int cindex;
   ec_slavet *csl;

   csl = &(context->slavelist[slave]);
   cindex = ec_findconfig( csl->eep_man, csl->eep_id );
   csl->configindex= cindex;
   /* slave found in configuration table ? */
   if (cindex)
   {
      csl->Dtype = ec_configlist[cindex].Dtype;
      strcpy(csl->name ,ec_configlist[cindex].name);
      csl->Ibits = ec_configlist[cindex].Ibits;
      csl->Obits = ec_configlist[cindex].Obits;
      if (csl->Obits)
      {
         csl->FMMU0func = 1;
      }
      if (csl->Ibits)
      {
         csl->FMMU1func = 2;
      }
      csl->FMMU[0].FMMUactive = ec_configlist[cindex].FM0ac;
      csl->FMMU[1].FMMUactive = ec_configlist[cindex].FM1ac;
      csl->SM[2].StartAddr = htoes(ec_configlist[cindex].SM2a);
      csl->SM[2].SMflags = htoel(ec_configlist[cindex].SM2f);
      /* simple (no mailbox) output slave found ? */
      if (csl->Obits && !csl->SM[2].StartAddr)
      {
         csl->SM[0].StartAddr = htoes(0x0f00);
         csl->SM[0].SMlength = htoes((csl->Obits + 7) / 8);
         csl->SM[0].SMflags = htoel(EC_DEFAULTDOSM0);
         csl->FMMU[0].FMMUactive = 1;
         csl->FMMU[0].FMMUtype = 2;
         csl->SMtype[0] = 3;
      }
      /* complex output slave */
      else
      {
         csl->SM[2].SMlength = htoes((csl->Obits + 7) / 8);
         csl->SMtype[2] = 3;
      }
      csl->SM[3].StartAddr = htoes(ec_configlist[cindex].SM3a);
      csl->SM[3].SMflags = htoel(ec_configlist[cindex].SM3f);
      /* simple (no mailbox) input slave found ? */
      if (csl->Ibits && !csl->SM[3].StartAddr)
      {
         csl->SM[1].StartAddr = htoes(0x1000);
         csl->SM[1].SMlength = htoes((csl->Ibits + 7) / 8);
         csl->SM[1].SMflags = htoel(0x00000000);
         csl->FMMU[1].FMMUactive = 1;
         csl->FMMU[1].FMMUtype = 1;
         csl->SMtype[1] = 4;
      }
      /* complex input slave */
      else
      {
         csl->SM[3].SMlength = htoes((csl->Ibits + 7) / 8);
         csl->SMtype[3] = 4;
      }
   }
   return cindex;
}
#else
static int ecx_config_from_table(ecx_contextt *context, uint16 slave)
{
   return 0;
}
#endif

/* If slave has SII and same slave ID done before, use previous data.
 * This is safe because SII is constant for same slave ID.
 */
static int ecx_lookup_prev_sii(ecx_contextt *context, uint16 slave)
{
   int i, nSM;
   if ((slave > 1) && (*(context->slavecount) > 0))
   {
      i = 1;
      while(((context->slavelist[i].eep_man != context->slavelist[slave].eep_man) ||
             (context->slavelist[i].eep_id  != context->slavelist[slave].eep_id ) ||
             (context->slavelist[i].eep_rev != context->slavelist[slave].eep_rev)) &&
            (i < slave))
      {
         i++;
      }
      if(i < slave)
      {
         context->slavelist[slave].CoEdetails = context->slavelist[i].CoEdetails;
         context->slavelist[slave].FoEdetails = context->slavelist[i].FoEdetails;
         context->slavelist[slave].EoEdetails = context->slavelist[i].EoEdetails;
         context->slavelist[slave].SoEdetails = context->slavelist[i].SoEdetails;
         if(context->slavelist[i].blockLRW > 0)
         {
            context->slavelist[slave].blockLRW = 1;
            context->slavelist[0].blockLRW++;
         }
         context->slavelist[slave].Ebuscurrent = context->slavelist[i].Ebuscurrent;
         context->slavelist[0].Ebuscurrent += context->slavelist[slave].Ebuscurrent;
         memcpy(context->slavelist[slave].name, context->slavelist[i].name, EC_MAXNAME + 1);
         for( nSM=0 ; nSM < EC_MAXSM ; nSM++ )
         {
            context->slavelist[slave].SM[nSM].StartAddr = context->slavelist[i].SM[nSM].StartAddr;
            context->slavelist[slave].SM[nSM].SMlength  = context->slavelist[i].SM[nSM].SMlength;
            context->slavelist[slave].SM[nSM].SMflags   = context->slavelist[i].SM[nSM].SMflags;
         }
         context->slavelist[slave].FMMU0func = context->slavelist[i].FMMU0func;
         context->slavelist[slave].FMMU1func = context->slavelist[i].FMMU1func;
         context->slavelist[slave].FMMU2func = context->slavelist[i].FMMU2func;
         context->slavelist[slave].FMMU3func = context->slavelist[i].FMMU3func;
         EC_PRINT("Copy SII slave %d from %d.\n", slave, i);
         return 1;
      }
   }
   return 0;
}

/** Enumerate and init all slaves.
 *
 * @param[in] context      = context struct
 * @param[in] usetable     = TRUE when using configtable to init slaves, FALSE otherwise
 * @return Workcounter of slave discover datagram = number of slaves found
 */
int ecx_config_init(ecx_contextt *context, uint8 usetable)
{
   uint16 slave, ADPh, configadr, ssigen;
   uint16 topology, estat;
   int16 topoc, slavec, aliasadr;
   uint8 b,h;
   uint8 SMc;
   uint32 eedat;
   int wkc, cindex, nSM;

   EC_PRINT("ec_config_init %d\n",usetable);
   ecx_init_context(context);
   wkc = ecx_detect_slaves(context);
   if (wkc > 0)
   {
      ecx_set_slaves_to_default(context);
      for (slave = 1; slave <= *(context->slavecount); slave++)
      {
         ADPh = (uint16)(1 - slave);
         context->slavelist[slave].Itype =
            etohs(ecx_APRDw(context->port, ADPh, ECT_REG_PDICTL, EC_TIMEOUTRET3)); /* read interface type of slave */
         /* a node offset is used to improve readibility of network frames */
         /* this has no impact on the number of addressable slaves (auto wrap around) */
         ecx_APWRw(context->port, ADPh, ECT_REG_STADR, htoes(slave + EC_NODEOFFSET) , EC_TIMEOUTRET3); /* set node address of slave */
         if (slave == 1)
         {
            b = 1; /* kill non ecat frames for first slave */
         }
         else
         {
            b = 0; /* pass all frames for following slaves */
         }
         ecx_APWRw(context->port, ADPh, ECT_REG_DLCTL, htoes(b), EC_TIMEOUTRET3); /* set non ecat frame behaviour */
         configadr = etohs(ecx_APRDw(context->port, ADPh, ECT_REG_STADR, EC_TIMEOUTRET3));
         context->slavelist[slave].configadr = configadr;
         ecx_FPRD(context->port, configadr, ECT_REG_ALIAS, sizeof(aliasadr), &aliasadr, EC_TIMEOUTRET3);
         context->slavelist[slave].aliasadr = etohs(aliasadr);
         ecx_FPRD(context->port, configadr, ECT_REG_EEPSTAT, sizeof(estat), &estat, EC_TIMEOUTRET3);
         estat = etohs(estat);
         if (estat & EC_ESTAT_R64) /* check if slave can read 8 byte chunks */
         {
            context->slavelist[slave].eep_8byte = 1;
         }
         ecx_readeeprom1(context, slave, ECT_SII_MANUF); /* Manuf */
      }
      for (slave = 1; slave <= *(context->slavecount); slave++)
      {
         context->slavelist[slave].eep_man =
            etohl(ecx_readeeprom2(context, slave, EC_TIMEOUTEEP)); /* Manuf */
         ecx_readeeprom1(context, slave, ECT_SII_ID); /* ID */
      }
      for (slave = 1; slave <= *(context->slavecount); slave++)
      {
         context->slavelist[slave].eep_id =
            etohl(ecx_readeeprom2(context, slave, EC_TIMEOUTEEP)); /* ID */
         ecx_readeeprom1(context, slave, ECT_SII_REV); /* revision */
      }
      for (slave = 1; slave <= *(context->slavecount); slave++)
      {
         context->slavelist[slave].eep_rev =
            etohl(ecx_readeeprom2(context, slave, EC_TIMEOUTEEP)); /* revision */
         ecx_readeeprom1(context, slave, ECT_SII_RXMBXADR); /* write mailbox address + mailboxsize */
      }
      for (slave = 1; slave <= *(context->slavecount); slave++)
      {
         eedat = etohl(ecx_readeeprom2(context, slave, EC_TIMEOUTEEP)); /* write mailbox address and mailboxsize */
         context->slavelist[slave].mbx_wo = (uint16)LO_WORD(eedat);
         context->slavelist[slave].mbx_l = (uint16)HI_WORD(eedat);
         if (context->slavelist[slave].mbx_l > 0)
         {
            ecx_readeeprom1(context, slave, ECT_SII_TXMBXADR); /* read mailbox offset */
         }
      }
      for (slave = 1; slave <= *(context->slavecount); slave++)
      {
         if (context->slavelist[slave].mbx_l > 0)
         {
            eedat = etohl(ecx_readeeprom2(context, slave, EC_TIMEOUTEEP)); /* read mailbox offset */
            context->slavelist[slave].mbx_ro = (uint16)LO_WORD(eedat); /* read mailbox offset */
            context->slavelist[slave].mbx_rl = (uint16)HI_WORD(eedat); /*read mailbox length */
            if (context->slavelist[slave].mbx_rl == 0)
            {
               context->slavelist[slave].mbx_rl = context->slavelist[slave].mbx_l;
            }
            ecx_readeeprom1(context, slave, ECT_SII_MBXPROTO);
         }
         configadr = context->slavelist[slave].configadr;
         if ((etohs(ecx_FPRDw(context->port, configadr, ECT_REG_ESCSUP, EC_TIMEOUTRET3)) & 0x04) > 0)  /* Support DC? */
         {
            context->slavelist[slave].hasdc = TRUE;
         }
         else
         {
            context->slavelist[slave].hasdc = FALSE;
         }
         topology = etohs(ecx_FPRDw(context->port, configadr, ECT_REG_DLSTAT, EC_TIMEOUTRET3)); /* extract topology from DL status */
         h = 0;
         b = 0;
         if ((topology & 0x0300) == 0x0200) /* port0 open and communication established */
         {
            h++;
            b |= 0x01;
         }
         if ((topology & 0x0c00) == 0x0800) /* port1 open and communication established */
         {
            h++;
            b |= 0x02;
         }
         if ((topology & 0x3000) == 0x2000) /* port2 open and communication established */
         {
            h++;
            b |= 0x04;
         }
         if ((topology & 0xc000) == 0x8000) /* port3 open and communication established */
         {
            h++;
            b |= 0x08;
         }
         /* ptype = Physical type*/
         context->slavelist[slave].ptype =
            LO_BYTE(etohs(ecx_FPRDw(context->port, configadr, ECT_REG_PORTDES, EC_TIMEOUTRET3)));
         context->slavelist[slave].topology = h;
         context->slavelist[slave].activeports = b;
         /* 0=no links, not possible             */
         /* 1=1 link  , end of line              */
         /* 2=2 links , one before and one after */
         /* 3=3 links , split point              */
         /* 4=4 links , cross point              */
         /* search for parent */
         context->slavelist[slave].parent = 0; /* parent is master */
         if (slave > 1)
         {
            topoc = 0;
            slavec = slave - 1;
            do
            {
               topology = context->slavelist[slavec].topology;
               if (topology == 1)
               {
                  topoc--; /* endpoint found */
               }
               if (topology == 3)
               {
                  topoc++; /* split found */
               }
               if (topology == 4)
               {
                  topoc += 2; /* cross found */
               }
               if (((topoc >= 0) && (topology > 1)) ||
                   (slavec == 1)) /* parent found */
               {
                  context->slavelist[slave].parent = slavec;
                  slavec = 1;
               }
               slavec--;
            }
            while (slavec > 0);
         }
         (void)ecx_statecheck(context, slave, EC_STATE_INIT,  EC_TIMEOUTSTATE); //* check state change Init */

         /* set default mailbox configuration if slave has mailbox */
         if (context->slavelist[slave].mbx_l>0)
         {
            context->slavelist[slave].SMtype[0] = 1;
            context->slavelist[slave].SMtype[1] = 2;
            context->slavelist[slave].SMtype[2] = 3;
            context->slavelist[slave].SMtype[3] = 4;
            context->slavelist[slave].SM[0].StartAddr = htoes(context->slavelist[slave].mbx_wo);
            context->slavelist[slave].SM[0].SMlength = htoes(context->slavelist[slave].mbx_l);
            context->slavelist[slave].SM[0].SMflags = htoel(EC_DEFAULTMBXSM0);
            context->slavelist[slave].SM[1].StartAddr = htoes(context->slavelist[slave].mbx_ro);
            context->slavelist[slave].SM[1].SMlength = htoes(context->slavelist[slave].mbx_rl);
            context->slavelist[slave].SM[1].SMflags = htoel(EC_DEFAULTMBXSM1);
            context->slavelist[slave].mbx_proto =
               ecx_readeeprom2(context, slave, EC_TIMEOUTEEP);
         }
         cindex = 0;
         /* use configuration table ? */
         if (usetable == 1)
         {
            cindex = ecx_config_from_table(context, slave);
         }
         /* slave not in configuration table, find out via SII */
         if (!cindex && !ecx_lookup_prev_sii(context, slave))
         {
            ssigen = ecx_siifind(context, slave, ECT_SII_GENERAL);
            /* SII general section */
            if (ssigen)
            {
               context->slavelist[slave].CoEdetails = ecx_siigetbyte(context, slave, ssigen + 0x07);
               context->slavelist[slave].FoEdetails = ecx_siigetbyte(context, slave, ssigen + 0x08);
               context->slavelist[slave].EoEdetails = ecx_siigetbyte(context, slave, ssigen + 0x09);
               context->slavelist[slave].SoEdetails = ecx_siigetbyte(context, slave, ssigen + 0x0a);
               if(1==1 || (ecx_siigetbyte(context, slave, ssigen + 0x0d) & 0x02) > 0)
               {
                  context->slavelist[slave].blockLRW = 1;
                  context->slavelist[0].blockLRW++;
               }
               context->slavelist[slave].Ebuscurrent = ecx_siigetbyte(context, slave, ssigen + 0x0e);
               context->slavelist[slave].Ebuscurrent += ecx_siigetbyte(context, slave, ssigen + 0x0f) << 8;
               context->slavelist[0].Ebuscurrent += context->slavelist[slave].Ebuscurrent;
            }
            /* SII strings section */
            if (ecx_siifind(context, slave, ECT_SII_STRING) > 0)
            {
               ecx_siistring(context, context->slavelist[slave].name, slave, 1);
            }
            /* no name for slave found, use constructed name */
            else
            {
               sprintf(context->slavelist[slave].name, "? M:%8.8x I:%8.8x",
                       (unsigned int)context->slavelist[slave].eep_man,
                       (unsigned int)context->slavelist[slave].eep_id);
            }
            /* SII SM section */
            nSM = ecx_siiSM(context, slave, context->eepSM);
            if (nSM>0)
            {
               context->slavelist[slave].SM[0].StartAddr = htoes(context->eepSM->PhStart);
               context->slavelist[slave].SM[0].SMlength = htoes(context->eepSM->Plength);
               context->slavelist[slave].SM[0].SMflags =
                  htoel((context->eepSM->Creg) + (context->eepSM->Activate << 16));
               SMc = 1;
               while ((SMc < EC_MAXSM) &&  ecx_siiSMnext(context, slave, context->eepSM, SMc))
               {
                  context->slavelist[slave].SM[SMc].StartAddr = htoes(context->eepSM->PhStart);
                  context->slavelist[slave].SM[SMc].SMlength = htoes(context->eepSM->Plength);
                  context->slavelist[slave].SM[SMc].SMflags =
                     htoel((context->eepSM->Creg) + (context->eepSM->Activate << 16));
                  SMc++;
               }
            }
            /* SII FMMU section */
            if (ecx_siiFMMU(context, slave, context->eepFMMU))
            {
               if (context->eepFMMU->FMMU0 !=0xff)
               {
                  context->slavelist[slave].FMMU0func = context->eepFMMU->FMMU0;
               }
               if (context->eepFMMU->FMMU1 !=0xff)
               {
                  context->slavelist[slave].FMMU1func = context->eepFMMU->FMMU1;
               }
               if (context->eepFMMU->FMMU2 !=0xff)
               {
                  context->slavelist[slave].FMMU2func = context->eepFMMU->FMMU2;
               }
               if (context->eepFMMU->FMMU3 !=0xff)
               {
                  context->slavelist[slave].FMMU3func = context->eepFMMU->FMMU3;
               }
            }
         }

         if (context->slavelist[slave].mbx_l > 0)
         {
            if (context->slavelist[slave].SM[0].StartAddr == 0x0000) /* should never happen */
            {
               EC_PRINT("Slave %d has no proper mailbox in configuration, try default.\n", slave);
               context->slavelist[slave].SM[0].StartAddr = htoes(0x1000);
               context->slavelist[slave].SM[0].SMlength = htoes(0x0080);
               context->slavelist[slave].SM[0].SMflags = htoel(EC_DEFAULTMBXSM0);
               context->slavelist[slave].SMtype[0] = 1;
            }
            if (context->slavelist[slave].SM[1].StartAddr == 0x0000) /* should never happen */
            {
               EC_PRINT("Slave %d has no proper mailbox out configuration, try default.\n", slave);
               context->slavelist[slave].SM[1].StartAddr = htoes(0x1080);
               context->slavelist[slave].SM[1].SMlength = htoes(0x0080);
               context->slavelist[slave].SM[1].SMflags = htoel(EC_DEFAULTMBXSM1);
               context->slavelist[slave].SMtype[1] = 2;
            }
            /* program SM0 mailbox in and SM1 mailbox out for slave */
            /* writing both SM in one datagram will solve timing issue in old NETX */
            ecx_FPWR(context->port, configadr, ECT_REG_SM0, sizeof(ec_smt) * 2,
               &(context->slavelist[slave].SM[0]), EC_TIMEOUTRET3);
         }
         /* some slaves need eeprom available to PDI in init->preop transition */
         ecx_eeprom2pdi(context, slave);
         /* request pre_op for slave */
         ecx_FPWRw(context->port, configadr, ECT_REG_ALCTL, htoes(EC_STATE_PRE_OP | EC_STATE_ACK) , EC_TIMEOUTRET3); /* set preop status */
      }
   }
   return wkc;
}

/* If slave has SII mapping and same slave ID done before, use previous mapping.
 * This is safe because SII mapping is constant for same slave ID.
 */
static int ecx_lookup_mapping(ecx_contextt *context, uint16 slave, int *Osize, int *Isize)
{
   int i, nSM;
   if ((slave > 1) && (*(context->slavecount) > 0))
   {
      i = 1;
      while(((context->slavelist[i].eep_man != context->slavelist[slave].eep_man) ||
             (context->slavelist[i].eep_id  != context->slavelist[slave].eep_id ) ||
             (context->slavelist[i].eep_rev != context->slavelist[slave].eep_rev)) &&
            (i < slave))
      {
         i++;
      }
      if(i < slave)
      {
         for( nSM=0 ; nSM < EC_MAXSM ; nSM++ )
         {
            context->slavelist[slave].SM[nSM].SMlength = context->slavelist[i].SM[nSM].SMlength;
            context->slavelist[slave].SMtype[nSM] = context->slavelist[i].SMtype[nSM];
         }
         *Osize = context->slavelist[i].Obits;
         *Isize = context->slavelist[i].Ibits;
         context->slavelist[slave].Obits = *Osize;
         context->slavelist[slave].Ibits = *Isize;
         EC_PRINT("Copy mapping slave %d from %d.\n", slave, i);
         return 1;
      }
   }
   return 0;
}

static int ecx_map_coe_soe(ecx_contextt *context, uint16 slave, int thread_n)
{
   int Isize, Osize;
   int rval;

   ecx_statecheck(context, slave, EC_STATE_PRE_OP, EC_TIMEOUTSTATE); /* check state change pre-op */

   EC_PRINT(" >Slave %d, configadr %x, state %2.2x\n",
            slave, context->slavelist[slave].configadr, context->slavelist[slave].state);

   /* execute special slave configuration hook Pre-Op to Safe-OP */
   if(context->slavelist[slave].PO2SOconfig) /* only if registered */
   {
      context->slavelist[slave].PO2SOconfig(slave);
   }
   /* if slave not found in configlist find IO mapping in slave self */
   if (!context->slavelist[slave].configindex)
   {
      Isize = 0;
      Osize = 0;
      if (context->slavelist[slave].mbx_proto & ECT_MBXPROT_COE) /* has CoE */
      {
         rval = 0;
         if (context->slavelist[slave].CoEdetails & ECT_COEDET_SDOCA) /* has Complete Access */
         {
            /* read PDO mapping via CoE and use Complete Access */
            rval = ecx_readPDOmapCA(context, slave, thread_n, &Osize, &Isize);
         }
         if (!rval) /* CA not available or not succeeded */
         {
            /* read PDO mapping via CoE */
            rval = ecx_readPDOmap(context, slave, &Osize, &Isize);
         }
         EC_PRINT("  CoE Osize:%d Isize:%d\n", Osize, Isize);
      }
      if ((!Isize && !Osize) && (context->slavelist[slave].mbx_proto & ECT_MBXPROT_SOE)) /* has SoE */
      {
         /* read AT / MDT mapping via SoE */
         rval = ecx_readIDNmap(context, slave, &Osize, &Isize);
         context->slavelist[slave].SM[2].SMlength = htoes((Osize + 7) / 8);
         context->slavelist[slave].SM[3].SMlength = htoes((Isize + 7) / 8);
         EC_PRINT("  SoE Osize:%d Isize:%d\n", Osize, Isize);
      }
      context->slavelist[slave].Obits = Osize;
      context->slavelist[slave].Ibits = Isize;
   }

   return 1;
}

static int ecx_map_sii(ecx_contextt *context, uint16 slave)
{
   int Isize, Osize;
   int nSM;
   ec_eepromPDOt eepPDO;

   Osize = context->slavelist[slave].Obits;
   Isize = context->slavelist[slave].Ibits;

   if (!Isize && !Osize) /* find PDO in previous slave with same ID */
   {
      (void)ecx_lookup_mapping(context, slave, &Osize, &Isize);
   }
   if (!Isize && !Osize) /* find PDO mapping by SII */
   {
      memset(&eepPDO, 0, sizeof(eepPDO));
      Isize = (int)ecx_siiPDO(context, slave, &eepPDO, 0);
      EC_PRINT("  SII Isize:%d\n", Isize);
      for( nSM=0 ; nSM < EC_MAXSM ; nSM++ )
      {
         if (eepPDO.SMbitsize[nSM] > 0)
         {
            context->slavelist[slave].SM[nSM].SMlength =  htoes((eepPDO.SMbitsize[nSM] + 7) / 8);
            context->slavelist[slave].SMtype[nSM] = 4;
            EC_PRINT("    SM%d length %d\n", nSM, eepPDO.SMbitsize[nSM]);
         }
      }
      Osize = (int)ecx_siiPDO(context, slave, &eepPDO, 1);
      EC_PRINT("  SII Osize:%d\n", Osize);
      for( nSM=0 ; nSM < EC_MAXSM ; nSM++ )
      {
         if (eepPDO.SMbitsize[nSM] > 0)
         {
            context->slavelist[slave].SM[nSM].SMlength =  htoes((eepPDO.SMbitsize[nSM] + 7) / 8);
            context->slavelist[slave].SMtype[nSM] = 3;
            EC_PRINT("    SM%d length %d\n", nSM, eepPDO.SMbitsize[nSM]);
         }
      }
   }
   context->slavelist[slave].Obits = Osize;
   context->slavelist[slave].Ibits = Isize;
   EC_PRINT("     ISIZE:%d %d OSIZE:%d\n",
      context->slavelist[slave].Ibits, Isize,context->slavelist[slave].Obits);

   return 1;
}

static int ecx_map_sm(ecx_contextt *context, uint16 slave)
{
   uint16 configadr;
   int nSM;

   configadr = context->slavelist[slave].configadr;

   EC_PRINT("  SM programming\n");
   if (!context->slavelist[slave].mbx_l && context->slavelist[slave].SM[0].StartAddr)
   {
      ecx_FPWR(context->port, configadr, ECT_REG_SM0,
         sizeof(ec_smt), &(context->slavelist[slave].SM[0]), EC_TIMEOUTRET3);
      EC_PRINT("    SM0 Type:%d StartAddr:%4.4x Flags:%8.8x\n",
          context->slavelist[slave].SMtype[0],
          context->slavelist[slave].SM[0].StartAddr,
          context->slavelist[slave].SM[0].SMflags);
   }
   if (!context->slavelist[slave].mbx_l && context->slavelist[slave].SM[1].StartAddr)
   {
      ecx_FPWR(context->port, configadr, ECT_REG_SM1,
         sizeof(ec_smt), &context->slavelist[slave].SM[1], EC_TIMEOUTRET3);
      EC_PRINT("    SM1 Type:%d StartAddr:%4.4x Flags:%8.8x\n",
          context->slavelist[slave].SMtype[1],
          context->slavelist[slave].SM[1].StartAddr,
          context->slavelist[slave].SM[1].SMflags);
   }
   /* program SM2 to SMx */
   for( nSM = 2 ; nSM < EC_MAXSM ; nSM++ )
   {
      if (context->slavelist[slave].SM[nSM].StartAddr)
      {
         /* check if SM length is zero -> clear enable flag */
         if( context->slavelist[slave].SM[nSM].SMlength == 0)
         {
            context->slavelist[slave].SM[nSM].SMflags =
               htoel( etohl(context->slavelist[slave].SM[nSM].SMflags) & EC_SMENABLEMASK);
         }
         ecx_FPWR(context->port, configadr, ECT_REG_SM0 + (nSM * sizeof(ec_smt)),
            sizeof(ec_smt), &context->slavelist[slave].SM[nSM], EC_TIMEOUTRET3);
         EC_PRINT("    SM%d Type:%d StartAddr:%4.4x Flags:%8.8x\n", nSM,
             context->slavelist[slave].SMtype[nSM],
             context->slavelist[slave].SM[nSM].StartAddr,
             context->slavelist[slave].SM[nSM].SMflags);
      }
   }
   if (context->slavelist[slave].Ibits > 7)
   {
      context->slavelist[slave].Ibytes = (context->slavelist[slave].Ibits + 7) / 8;
   }
   if (context->slavelist[slave].Obits > 7)
   {
      context->slavelist[slave].Obytes = (context->slavelist[slave].Obits + 7) / 8;
   }

   return 1;
}

OSAL_THREAD_FUNC ecx_mapper_thread(void *param)
{
   ecx_mapt_t *maptp;
   maptp = param;
   ecx_map_coe_soe(maptp->context, maptp->slave, maptp->thread_n);
   maptp->running = 0;
}

static int ecx_find_mapt(void)
{
   int p;
   p = 0;
   while((p < EC_MAX_MAPT) && ecx_mapt[p].running)
   {
      p++;
   }
   if(p < EC_MAX_MAPT)
   {
      return p;
   }
   else
   {
      return -1;
   }
}

static int ecx_get_threadcount(void)
{
   int thrc, thrn;
   thrc = 0;
   for(thrn = 0 ; thrn < EC_MAX_MAPT ; thrn++)
   {
      thrc += ecx_mapt[thrn].running;
   }
   return thrc;
}

/** Map all PDOs in one group of slaves to IOmap.
 *
 * @param[in]  context    = context struct
 * @param[out] pIOmap     = pointer to IOmap
 * @param[in]  group      = group to map, 0 = all groups
 * @return IOmap size
 */
int ecx_config_map_group(ecx_contextt *context, void *pIOmap, uint8 group)
{
   uint16 slave, configadr;
   int BitCount, ByteCount, FMMUsize, FMMUdone;
   uint16 SMlength, EndAddr;
   uint8 BitPos;
   uint8 SMc, FMMUc;
   uint32 LogAddr = 0;
   uint32 oLogAddr = 0;
   uint32 diff;
   uint16 currentsegment = 0;
   uint32 segmentsize = 0;
   int thrn, thrc;

   if ((*(context->slavecount) > 0) && (group < context->maxgroup))
   {
      EC_PRINT("ec_config_map_group IOmap:%p group:%d\n", pIOmap, group);
      LogAddr = context->grouplist[group].logstartaddr;
      oLogAddr = LogAddr;
      BitPos = 0;
      context->grouplist[group].nsegments = 0;
      context->grouplist[group].outputsWKC = 0;
      context->grouplist[group].inputsWKC = 0;

      for(thrn = 0 ; thrn < EC_MAX_MAPT ; thrn++)
      {
         ecx_mapt[thrn].running = 0;
      }
      /* find CoE and SoE mapping of slaves in multiple threads */
      for (slave = 1; slave <= *(context->slavecount); slave++)
      {
         if (!group || (group == context->slavelist[slave].group))
         {
            if(EC_MAX_MAPT <= 1)
            {
               /* serialised version */
               ecx_map_coe_soe(context, slave, 0);
            }
            else
            {
               /* multi-threaded version */
               while((thrn = ecx_find_mapt()) < 0)
               {
                  osal_usleep(1000);
               }
               ecx_mapt[thrn].context = context;
               ecx_mapt[thrn].slave = slave;
               ecx_mapt[thrn].thread_n = thrn;
               ecx_mapt[thrn].running = 1;
               osal_thread_create(&(ecx_threadh[thrn]), 128000,
                  &ecx_mapper_thread, &(ecx_mapt[thrn]));
            }
         }
      }
      /* wait for all threads to finish */
      do
      {
         thrc = ecx_get_threadcount();
         if(thrc)
         {
            osal_usleep(1000);
         }
      } while(thrc);
      /* find SII mapping of slave and program SM */
      for (slave = 1; slave <= *(context->slavecount); slave++)
      {
         if (!group || (group == context->slavelist[slave].group))
         {
            ecx_map_sii(context, slave);
            ecx_map_sm(context, slave);
         }
      }

      /* do input mapping of slave and program FMMUs */
      for (slave = 1; slave <= *(context->slavecount); slave++)
      {
         configadr = context->slavelist[slave].configadr;

         if (!group || (group == context->slavelist[slave].group))
         {
            FMMUc = context->slavelist[slave].FMMUunused;
            SMc = 0;
            BitCount = 0;
            ByteCount = 0;
            EndAddr = 0;
            FMMUsize = 0;
            FMMUdone = 0;
            /* create output mapping */
            if (context->slavelist[slave].Obits)
            {
               EC_PRINT("  OUTPUT MAPPING\n");
               /* search for SM that contribute to the output mapping */
               while ( (SMc < (EC_MAXSM - 1)) && (FMMUdone < ((context->slavelist[slave].Obits + 7) / 8)))
               {
                  EC_PRINT("    FMMU %d\n", FMMUc);
                  while ( (SMc < (EC_MAXSM - 1)) && (context->slavelist[slave].SMtype[SMc] != 3)) SMc++;
                  EC_PRINT("      SM%d\n", SMc);
                  context->slavelist[slave].FMMU[FMMUc].PhysStart =
                     context->slavelist[slave].SM[SMc].StartAddr;
                  SMlength = etohs(context->slavelist[slave].SM[SMc].SMlength);
                  ByteCount += SMlength;
                  BitCount += SMlength * 8;
                  EndAddr = etohs(context->slavelist[slave].SM[SMc].StartAddr) + SMlength;
                  while ( (BitCount < context->slavelist[slave].Obits) && (SMc < (EC_MAXSM - 1)) ) /* more SM for output */
                  {
                     SMc++;
                     while ( (SMc < (EC_MAXSM - 1)) && (context->slavelist[slave].SMtype[SMc] != 3)) SMc++;
                     /* if addresses from more SM connect use one FMMU otherwise break up in mutiple FMMU */
                     if ( etohs(context->slavelist[slave].SM[SMc].StartAddr) > EndAddr )
                     {
                        break;
                     }
                     EC_PRINT("      SM%d\n", SMc);
                     SMlength = etohs(context->slavelist[slave].SM[SMc].SMlength);
                     ByteCount += SMlength;
                     BitCount += SMlength * 8;
                     EndAddr = etohs(context->slavelist[slave].SM[SMc].StartAddr) + SMlength;
                  }

                  /* bit oriented slave */
                  if (!context->slavelist[slave].Obytes)
                  {
                     context->slavelist[slave].FMMU[FMMUc].LogStart = htoel(LogAddr);
                     context->slavelist[slave].FMMU[FMMUc].LogStartbit = BitPos;
                     BitPos += context->slavelist[slave].Obits - 1;
                     if (BitPos > 7)
                     {
                        LogAddr++;
                        BitPos -= 8;
                     }
                     FMMUsize = LogAddr - etohl(context->slavelist[slave].FMMU[FMMUc].LogStart) + 1;
                     context->slavelist[slave].FMMU[FMMUc].LogLength = htoes(FMMUsize);
                     context->slavelist[slave].FMMU[FMMUc].LogEndbit = BitPos;
                     BitPos ++;
                     if (BitPos > 7)
                     {
                        LogAddr++;
                        BitPos -= 8;
                     }
                  }
                  /* byte oriented slave */
                  else
                  {
                     if (BitPos)
                     {
                        LogAddr++;
                        BitPos = 0;
                     }
                     context->slavelist[slave].FMMU[FMMUc].LogStart = htoel(LogAddr);
                     context->slavelist[slave].FMMU[FMMUc].LogStartbit = BitPos;
                     BitPos = 7;
                     FMMUsize = ByteCount;
                     if ((FMMUsize + FMMUdone)> (int)context->slavelist[slave].Obytes)
                     {
                        FMMUsize = context->slavelist[slave].Obytes - FMMUdone;
                     }
                     LogAddr += FMMUsize;
                     context->slavelist[slave].FMMU[FMMUc].LogLength = htoes(FMMUsize);
                     context->slavelist[slave].FMMU[FMMUc].LogEndbit = BitPos;
                     BitPos = 0;
                  }
                  FMMUdone += FMMUsize;
                  context->slavelist[slave].FMMU[FMMUc].PhysStartBit = 0;
                  context->slavelist[slave].FMMU[FMMUc].FMMUtype = 2;
                  context->slavelist[slave].FMMU[FMMUc].FMMUactive = 1;
                  /* program FMMU for output */
                  ecx_FPWR(context->port, configadr, ECT_REG_FMMU0 + (sizeof(ec_fmmut) * FMMUc),
                     sizeof(ec_fmmut), &(context->slavelist[slave].FMMU[FMMUc]), EC_TIMEOUTRET3);
                  context->grouplist[group].outputsWKC++;
                  if (!context->slavelist[slave].outputs)
                  {
                     context->slavelist[slave].outputs =
                        (uint8 *)(pIOmap) + etohl(context->slavelist[slave].FMMU[FMMUc].LogStart);
                     context->slavelist[slave].Ostartbit =
                        context->slavelist[slave].FMMU[FMMUc].LogStartbit;
                     EC_PRINT("    slave %d Outputs %p startbit %d\n",
                        slave,
                        context->slavelist[slave].outputs,
                        context->slavelist[slave].Ostartbit);
                  }
                  FMMUc++;
               }
               context->slavelist[slave].FMMUunused = FMMUc;
               diff = LogAddr - oLogAddr;
               oLogAddr = LogAddr;
               if ((segmentsize + diff) > (EC_MAXLRWDATA - EC_FIRSTDCDATAGRAM))
               {
                  context->grouplist[group].IOsegment[currentsegment] = segmentsize;
                  if (currentsegment < (EC_MAXIOSEGMENTS - 1))
                  {
                     currentsegment++;
                     segmentsize = diff;
                  }
               }
               else
               {
                  segmentsize += diff;
               }
            }
         }
      }
      if (BitPos)
      {
         LogAddr++;
         oLogAddr = LogAddr;
         BitPos = 0;
         if ((segmentsize + 1) > (EC_MAXLRWDATA - EC_FIRSTDCDATAGRAM))
         {
            context->grouplist[group].IOsegment[currentsegment] = segmentsize;
            if (currentsegment < (EC_MAXIOSEGMENTS - 1))
            {
               currentsegment++;
               segmentsize = 1;
            }
         }
         else
         {
            segmentsize += 1;
         }
      }
      context->grouplist[group].outputs = pIOmap;
      context->grouplist[group].Obytes = LogAddr;
      context->grouplist[group].nsegments = currentsegment + 1;
      context->grouplist[group].Isegment = currentsegment;
      context->grouplist[group].Ioffset = segmentsize;
      if (!group)
      {
         context->slavelist[0].outputs = pIOmap;
         context->slavelist[0].Obytes = LogAddr; /* store output bytes in master record */
      }

      /* do input mapping of slave and program FMMUs */
      for (slave = 1; slave <= *(context->slavecount); slave++)
      {
         configadr = context->slavelist[slave].configadr;
         if (!group || (group == context->slavelist[slave].group))
         {
            FMMUc = context->slavelist[slave].FMMUunused;
            if (context->slavelist[slave].Obits) /* find free FMMU */
            {
               while ( context->slavelist[slave].FMMU[FMMUc].LogStart ) FMMUc++;
            }
            SMc = 0;
            BitCount = 0;
            ByteCount = 0;
            EndAddr = 0;
            FMMUsize = 0;
            FMMUdone = 0;
            /* create input mapping */
            if (context->slavelist[slave].Ibits)
            {
               EC_PRINT(" =Slave %d, INPUT MAPPING\n", slave);
               /* search for SM that contribute to the input mapping */
               while ( (SMc < (EC_MAXSM - 1)) && (FMMUdone < ((context->slavelist[slave].Ibits + 7) / 8)))
               {
                  EC_PRINT("    FMMU %d\n", FMMUc);
                  while ( (SMc < (EC_MAXSM - 1)) && (context->slavelist[slave].SMtype[SMc] != 4)) SMc++;
                  EC_PRINT("      SM%d\n", SMc);
                  context->slavelist[slave].FMMU[FMMUc].PhysStart =
                     context->slavelist[slave].SM[SMc].StartAddr;
                  SMlength = etohs(context->slavelist[slave].SM[SMc].SMlength);
                  ByteCount += SMlength;
                  BitCount += SMlength * 8;
                  EndAddr = etohs(context->slavelist[slave].SM[SMc].StartAddr) + SMlength;
                  while ( (BitCount < context->slavelist[slave].Ibits) && (SMc < (EC_MAXSM - 1)) ) /* more SM for input */
                  {
                     SMc++;
                     while ( (SMc < (EC_MAXSM - 1)) && (context->slavelist[slave].SMtype[SMc] != 4)) SMc++;
                     /* if addresses from more SM connect use one FMMU otherwise break up in mutiple FMMU */
                     if ( etohs(context->slavelist[slave].SM[SMc].StartAddr) > EndAddr )
                     {
                        break;
                     }
                     EC_PRINT("      SM%d\n", SMc);
                     SMlength = etohs(context->slavelist[slave].SM[SMc].SMlength);
                     ByteCount += SMlength;
                     BitCount += SMlength * 8;
                     EndAddr = etohs(context->slavelist[slave].SM[SMc].StartAddr) + SMlength;
                  }

                  /* bit oriented slave */
                  if (!context->slavelist[slave].Ibytes)
                  {
                     context->slavelist[slave].FMMU[FMMUc].LogStart = htoel(LogAddr);
                     context->slavelist[slave].FMMU[FMMUc].LogStartbit = BitPos;
                     BitPos += context->slavelist[slave].Ibits - 1;
                     if (BitPos > 7)
                     {
                        LogAddr++;
                        BitPos -= 8;
                     }
                     FMMUsize = LogAddr - etohl(context->slavelist[slave].FMMU[FMMUc].LogStart) + 1;
                     context->slavelist[slave].FMMU[FMMUc].LogLength = htoes(FMMUsize);
                     context->slavelist[slave].FMMU[FMMUc].LogEndbit = BitPos;
                     BitPos ++;
                     if (BitPos > 7)
                     {
                        LogAddr++;
                        BitPos -= 8;
                     }
                  }
                  /* byte oriented slave */
                  else
                  {
                     if (BitPos)
                     {
                        LogAddr++;
                        BitPos = 0;
                     }
                     context->slavelist[slave].FMMU[FMMUc].LogStart = htoel(LogAddr);
                     context->slavelist[slave].FMMU[FMMUc].LogStartbit = BitPos;
                     BitPos = 7;
                     FMMUsize = ByteCount;
                     if ((FMMUsize + FMMUdone)> (int)context->slavelist[slave].Ibytes)
                     {
                        FMMUsize = context->slavelist[slave].Ibytes - FMMUdone;
                     }
                     LogAddr += FMMUsize;
                     context->slavelist[slave].FMMU[FMMUc].LogLength = htoes(FMMUsize);
                     context->slavelist[slave].FMMU[FMMUc].LogEndbit = BitPos;
                     BitPos = 0;
                  }
                  FMMUdone += FMMUsize;
                  if (context->slavelist[slave].FMMU[FMMUc].LogLength)
                  {
                     context->slavelist[slave].FMMU[FMMUc].PhysStartBit = 0;
                     context->slavelist[slave].FMMU[FMMUc].FMMUtype = 1;
                     context->slavelist[slave].FMMU[FMMUc].FMMUactive = 1;
                     /* program FMMU for input */
                     ecx_FPWR(context->port, configadr, ECT_REG_FMMU0 + (sizeof(ec_fmmut) * FMMUc),
                        sizeof(ec_fmmut), &(context->slavelist[slave].FMMU[FMMUc]), EC_TIMEOUTRET3);
                     /* add one for an input FMMU */
                     context->grouplist[group].inputsWKC++;
                  }
                  if (!context->slavelist[slave].inputs)
                  {
                     context->slavelist[slave].inputs =
                        (uint8 *)(pIOmap) + etohl(context->slavelist[slave].FMMU[FMMUc].LogStart);
                     context->slavelist[slave].Istartbit =
                        context->slavelist[slave].FMMU[FMMUc].LogStartbit;
                     EC_PRINT("    Inputs %p startbit %d\n",
                        context->slavelist[slave].inputs,
                        context->slavelist[slave].Istartbit);
                  }
                  FMMUc++;
               }
               context->slavelist[slave].FMMUunused = FMMUc;
               diff = LogAddr - oLogAddr;
               oLogAddr = LogAddr;
               if ((segmentsize + diff) > (EC_MAXLRWDATA - EC_FIRSTDCDATAGRAM))
               {
                  context->grouplist[group].IOsegment[currentsegment] = segmentsize;
                  if (currentsegment < (EC_MAXIOSEGMENTS - 1))
                  {
                     currentsegment++;
                     segmentsize = diff;
                  }
               }
               else
               {
                  segmentsize += diff;
               }
            }

            ecx_eeprom2pdi(context, slave); /* set Eeprom control to PDI */
            ecx_FPWRw(context->port, configadr, ECT_REG_ALCTL, htoes(EC_STATE_SAFE_OP) , EC_TIMEOUTRET3); /* set safeop status */

            if (context->slavelist[slave].blockLRW)
            {
               context->grouplist[group].blockLRW++;
            }
            context->grouplist[group].Ebuscurrent += context->slavelist[slave].Ebuscurrent;
         }
      }
      if (BitPos)
      {
         LogAddr++;
         oLogAddr = LogAddr;
         BitPos = 0;
         if ((segmentsize + 1) > (EC_MAXLRWDATA - EC_FIRSTDCDATAGRAM))
         {
            context->grouplist[group].IOsegment[currentsegment] = segmentsize;
            if (currentsegment < (EC_MAXIOSEGMENTS - 1))
            {
               currentsegment++;
               segmentsize = 1;
            }
         }
         else
         {
            segmentsize += 1;
         }
      }
      context->grouplist[group].IOsegment[currentsegment] = segmentsize;
      context->grouplist[group].nsegments = currentsegment + 1;
      context->grouplist[group].inputs = (uint8 *)(pIOmap) + context->grouplist[group].Obytes;
      context->grouplist[group].Ibytes = LogAddr - context->grouplist[group].Obytes;
      if (!group)
      {
         context->slavelist[0].inputs = (uint8 *)(pIOmap) + context->slavelist[0].Obytes;
         context->slavelist[0].Ibytes = LogAddr - context->slavelist[0].Obytes; /* store input bytes in master record */
      }

      EC_PRINT("IOmapSize %d\n", LogAddr - context->grouplist[group].logstartaddr);

      return (LogAddr - context->grouplist[group].logstartaddr);
   }

   return 0;
}

/** Recover slave.
 *
 * @param[in] context = context struct
 * @param[in] slave   = slave to recover
 * @param[in] timeout = local timeout f.e. EC_TIMEOUTRET3
 * @return >0 if successful
 */
int ecx_recover_slave(ecx_contextt *context, uint16 slave, int timeout)
{
   int rval;
   int wkc;
   uint16 ADPh, configadr, readadr;

   rval = 0;
   configadr = context->slavelist[slave].configadr;
   ADPh = (uint16)(1 - slave);
   /* check if we found another slave than the requested */
   readadr = 0xfffe;
   wkc = ecx_APRD(context->port, ADPh, ECT_REG_STADR, sizeof(readadr), &readadr, timeout);
   /* correct slave found, finished */
   if(readadr == configadr)
   {
       return 1;
   }
   /* only try if no config address*/
   if( (wkc > 0) && (readadr == 0))
   {
      /* clear possible slaves at EC_TEMPNODE */
      ecx_FPWRw(context->port, EC_TEMPNODE, ECT_REG_STADR, htoes(0) , 0);
      /* set temporary node address of slave */
      if(ecx_APWRw(context->port, ADPh, ECT_REG_STADR, htoes(EC_TEMPNODE) , timeout) <= 0)
      {
         ecx_FPWRw(context->port, EC_TEMPNODE, ECT_REG_STADR, htoes(0) , 0);
         return 0; /* slave fails to respond */
      }

      context->slavelist[slave].configadr = EC_TEMPNODE; /* temporary config address */
      ecx_eeprom2master(context, slave); /* set Eeprom control to master */

      /* check if slave is the same as configured before */
      if ((ecx_FPRDw(context->port, EC_TEMPNODE, ECT_REG_ALIAS, timeout) ==
             context->slavelist[slave].aliasadr) &&
          (ecx_readeeprom(context, slave, ECT_SII_ID, EC_TIMEOUTEEP) ==
             context->slavelist[slave].eep_id) &&
          (ecx_readeeprom(context, slave, ECT_SII_MANUF, EC_TIMEOUTEEP) ==
             context->slavelist[slave].eep_man) &&
          (ecx_readeeprom(context, slave, ECT_SII_REV, EC_TIMEOUTEEP) ==
             context->slavelist[slave].eep_rev))
      {
         rval = ecx_FPWRw(context->port, EC_TEMPNODE, ECT_REG_STADR, htoes(configadr) , timeout);
         context->slavelist[slave].configadr = configadr;
      }
      else
      {
         /* slave is not the expected one, remove config address*/
         ecx_FPWRw(context->port, EC_TEMPNODE, ECT_REG_STADR, htoes(0) , timeout);
         context->slavelist[slave].configadr = configadr;
      }
   }

   return rval;
}

/** Reconfigure slave.
 *
 * @param[in] context = context struct
 * @param[in] slave   = slave to reconfigure
 * @param[in] timeout = local timeout f.e. EC_TIMEOUTRET3
 * @return Slave state
 */
int ecx_reconfig_slave(ecx_contextt *context, uint16 slave, int timeout)
{
   int state, nSM, FMMUc;
   uint16 configadr;

   configadr = context->slavelist[slave].configadr;
   if (ecx_FPWRw(context->port, configadr, ECT_REG_ALCTL, htoes(EC_STATE_INIT) , timeout) <= 0)
   {
      return 0;
   }
   state = 0;
   ecx_eeprom2pdi(context, slave); /* set Eeprom control to PDI */
   /* check state change init */
   state = ecx_statecheck(context, slave, EC_STATE_INIT, EC_TIMEOUTSTATE);
   if(state == EC_STATE_INIT)
   {
      /* program all enabled SM */
      for( nSM = 0 ; nSM < EC_MAXSM ; nSM++ )
      {
         if (context->slavelist[slave].SM[nSM].StartAddr)
         {
            ecx_FPWR(context->port, configadr, ECT_REG_SM0 + (nSM * sizeof(ec_smt)),
               sizeof(ec_smt), &context->slavelist[slave].SM[nSM], timeout);
         }
      }
      ecx_FPWRw(context->port, configadr, ECT_REG_ALCTL, htoes(EC_STATE_PRE_OP) , timeout);
      state = ecx_statecheck(context, slave, EC_STATE_PRE_OP, EC_TIMEOUTSTATE); /* check state change pre-op */
      if( state == EC_STATE_PRE_OP)
      {
         /* execute special slave configuration hook Pre-Op to Safe-OP */
         if(context->slavelist[slave].PO2SOconfig) /* only if registered */
         {
            context->slavelist[slave].PO2SOconfig(slave);
         }
         ecx_FPWRw(context->port, configadr, ECT_REG_ALCTL, htoes(EC_STATE_SAFE_OP) , timeout); /* set safeop status */
         state = ecx_statecheck(context, slave, EC_STATE_SAFE_OP, EC_TIMEOUTSTATE); /* check state change safe-op */
         /* program configured FMMU */
         for( FMMUc = 0 ; FMMUc < context->slavelist[slave].FMMUunused ; FMMUc++ )
         {
            ecx_FPWR(context->port, configadr, ECT_REG_FMMU0 + (sizeof(ec_fmmut) * FMMUc),
               sizeof(ec_fmmut), &context->slavelist[slave].FMMU[FMMUc], timeout);
         }
      }
   }

   return state;
}

#ifdef EC_VER1
/** Enumerate and init all slaves.
 *
 * @param[in] usetable     = TRUE when using configtable to init slaves, FALSE otherwise
 * @return Workcounter of slave discover datagram = number of slaves found
 * @see ecx_config_init
 */
int ec_config_init(uint8 usetable)
{
   return ecx_config_init(&ecx_context, usetable);
}

/** Map all PDOs in one group of slaves to IOmap.
 *
 * @param[out] pIOmap     = pointer to IOmap
 * @param[in]  group      = group to map, 0 = all groups
 * @return IOmap size
 * @see ecx_config_map_group
 */
int ec_config_map_group(void *pIOmap, uint8 group)
{
   return ecx_config_map_group(&ecx_context, pIOmap, group);
}

/** Map all PDOs from slaves to IOmap.
 *
 * @param[out] pIOmap     = pointer to IOmap
 * @return IOmap size
 */
int ec_config_map(void *pIOmap)
{
   return ec_config_map_group(pIOmap, 0);
}

/** Enumerate / map and init all slaves.
 *
 * @param[in] usetable    = TRUE when using configtable to init slaves, FALSE otherwise
 * @param[out] pIOmap     = pointer to IOmap
 * @return Workcounter of slave discover datagram = number of slaves found
 */
int ec_config(uint8 usetable, void *pIOmap)
{
   int wkc;
   wkc = ec_config_init(usetable);
   if (wkc)
   {
      ec_config_map(pIOmap);
   }
   return wkc;
}

/** Recover slave.
 *
 * @param[in] slave   = slave to recover
 * @param[in] timeout = local timeout f.e. EC_TIMEOUTRET3
 * @return >0 if successful
 * @see ecx_recover_slave
 */
int ec_recover_slave(uint16 slave, int timeout)
{
   return ecx_recover_slave(&ecx_context, slave, timeout);
}

/** Reconfigure slave.
 *
 * @param[in] slave   = slave to reconfigure
 * @param[in] timeout = local timeout f.e. EC_TIMEOUTRET3
 * @return Slave state
 * @see ecx_reconfig_slave
 */
int ec_reconfig_slave(uint16 slave, int timeout)
{
   return ecx_reconfig_slave(&ecx_context, slave, timeout);
}
#endif
