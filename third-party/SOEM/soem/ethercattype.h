/*
 * Simple Open EtherCAT Master Library
 *
 * File    : ethercattype.h
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
 * General typedefs and defines for EtherCAT.
 *
 * Defines that could need optimisation for specific applications
 * are the EC_TIMEOUTxxx. Assumptions for the standard settings are a
 * standard linux PC or laptop and a wired connection to maximal 100 slaves.
 * For use with wireless connections or lots of slaves the timeouts need
 * increasing. For fast systems running Xenomai and RT-net or alike the
 * timeouts need to be shorter.
 */

#ifndef _EC_TYPE_H
#define _EC_TYPE_H

#ifdef __cplusplus
extern "C"
{
#endif

/** Define Little or Big endian target */
#define EC_LITTLE_ENDIAN

/** define EC_VER1 if version 1 default context and functions are needed
 * comment if application uses only ecx_ functions and own context */
#define EC_VER1

#include "osal.h"

/** return value general error */
#define EC_ERROR           -3
/** return value no frame returned */
#define EC_NOFRAME         -1
/** return value unknown frame received */
#define EC_OTHERFRAME      -2
/** maximum EtherCAT frame length in bytes */
#define EC_MAXECATFRAME    1518
/** maximum EtherCAT LRW frame length in bytes */
/* MTU - Ethernet header - length - datagram header - WCK - FCS */
#define EC_MAXLRWDATA      (EC_MAXECATFRAME - 14 - 2 - 10 - 2 - 4)
/** size of DC datagram used in first LRW frame */
#define EC_FIRSTDCDATAGRAM 20
/** standard frame buffer size in bytes */
#define EC_BUFSIZE         EC_MAXECATFRAME
/** datagram type EtherCAT */
#define EC_ECATTYPE        0x1000
/** number of frame buffers per channel (tx, rx1 rx2) */
#define EC_MAXBUF          16
/** timeout value in us for tx frame to return to rx */
#define EC_TIMEOUTRET      2000
/** timeout value in us for safe data transfer, max. triple retry */
#define EC_TIMEOUTRET3     (EC_TIMEOUTRET * 3)
/** timeout value in us for return "safe" variant (f.e. wireless) */
#define EC_TIMEOUTSAFE     20000
/** timeout value in us for EEPROM access */
#define EC_TIMEOUTEEP      20000
/** timeout value in us for tx mailbox cycle */
#define EC_TIMEOUTTXM      20000
/** timeout value in us for rx mailbox cycle */
#define EC_TIMEOUTRXM      700000
/** timeout value in us for check statechange */
#define EC_TIMEOUTSTATE    2000000
/** size of EEPROM bitmap cache */
#define EC_MAXEEPBITMAP    128
/** size of EEPROM cache buffer */
#define EC_MAXEEPBUF       EC_MAXEEPBITMAP << 5
/** default number of retries if wkc <= 0 */
#define EC_DEFAULTRETRIES  3

/** definition for frame buffers */
typedef uint8 ec_bufT[EC_BUFSIZE];

/** ethernet header definition */
PACKED_BEGIN
typedef struct PACKED
{
   /** destination MAC */
   uint16  da0,da1,da2;
   /** source MAC */
   uint16  sa0,sa1,sa2;
   /** ethernet type */
   uint16  etype;
} ec_etherheadert;
PACKED_END

/** ethernet header size */
#define ETH_HEADERSIZE      sizeof(ec_etherheadert)

/** EtherCAT datagram header definition */
PACKED_BEGIN
typedef struct PACKED
{
   /** length of EtherCAT datagram */
   uint16  elength;
   /** EtherCAT command, see ec_cmdtype */
   uint8   command;
   /** index, used in SOEM for Tx to Rx recombination */
   uint8   index;
   /** ADP */
   uint16  ADP;
   /** ADO */
   uint16  ADO;
   /** length of data portion in datagram */
   uint16  dlength;
   /** interrupt, currently unused */
   uint16  irpt;
} ec_comt;
PACKED_END

/** EtherCAT header size */
#define EC_HEADERSIZE       sizeof(ec_comt)
/** size of ec_comt.elength item in EtherCAT header */
#define EC_ELENGTHSIZE      sizeof(uint16)
/** offset position of command in EtherCAT header */
#define EC_CMDOFFSET        EC_ELENGTHSIZE
/** size of workcounter item in EtherCAT datagram */
#define EC_WKCSIZE          sizeof(uint16)
/** definition of datagram follows bit in ec_comt.dlength */
#define EC_DATAGRAMFOLLOWS  (1 << 15)

/** Possible error codes returned. */
typedef enum
{
   /** No error */
   EC_ERR_OK         = 0,
   /** Library already initialized. */
   EC_ERR_ALREADY_INITIALIZED,
   /** Library not initialized. */
   EC_ERR_NOT_INITIALIZED,
   /** Timeout occurred during execution of the function. */
   EC_ERR_TIMEOUT,
   /** No slaves were found. */
   EC_ERR_NO_SLAVES,
   /** Function failed. */
   EC_ERR_NOK
} ec_err;

/** Possible EtherCAT slave states */
typedef enum
{
   /** Init state*/
   EC_STATE_INIT           = 0x01,
   /** Pre-operational. */
   EC_STATE_PRE_OP         = 0x02,
   /** Boot state*/
   EC_STATE_BOOT            = 0x03,
   /** Safe-operational. */
   EC_STATE_SAFE_OP        = 0x04,
   /** Operational */
   EC_STATE_OPERATIONAL    = 0x08,
   /** Error or ACK error */
   EC_STATE_ACK            = 0x10,
   EC_STATE_ERROR          = 0x10
} ec_state;

/** Possible buffer states */
typedef enum
{
   /** Empty */
   EC_BUF_EMPTY        = 0x00,
   /** Allocated, but not filled */
   EC_BUF_ALLOC        = 0x01,
   /** Transmitted */
   EC_BUF_TX           = 0x02,
   /** Received, but not consumed */
   EC_BUF_RCVD         = 0x03,
   /** Cycle completed */
   EC_BUF_COMPLETE     = 0x04
} ec_bufstate;

/** Ethercat data types */
typedef enum
{
   ECT_BOOLEAN         = 0x0001,
   ECT_INTEGER8        = 0x0002,
   ECT_INTEGER16       = 0x0003,
   ECT_INTEGER32       = 0x0004,
   ECT_UNSIGNED8       = 0x0005,
   ECT_UNSIGNED16      = 0x0006,
   ECT_UNSIGNED32      = 0x0007,
   ECT_REAL32          = 0x0008,
   ECT_VISIBLE_STRING  = 0x0009,
   ECT_OCTET_STRING    = 0x000A,
   ECT_UNICODE_STRING  = 0x000B,
   ECT_TIME_OF_DAY     = 0x000C,
   ECT_TIME_DIFFERENCE = 0x000D,
   ECT_DOMAIN          = 0x000F,
   ECT_INTEGER24       = 0x0010,
   ECT_REAL64          = 0x0011,
   ECT_INTEGER64       = 0x0015,
   ECT_UNSIGNED24      = 0x0016,
   ECT_UNSIGNED64      = 0x001B,
   ECT_BIT1            = 0x0030,
   ECT_BIT2            = 0x0031,
   ECT_BIT3            = 0x0032,
   ECT_BIT4            = 0x0033,
   ECT_BIT5            = 0x0034,
   ECT_BIT6            = 0x0035,
   ECT_BIT7            = 0x0036,
   ECT_BIT8            = 0x0037
} ec_datatype;

/** Ethercat command types */
typedef enum
{
   /** No operation */
   EC_CMD_NOP          = 0x00,
   /** Auto Increment Read */
   EC_CMD_APRD,
   /** Auto Increment Write */
   EC_CMD_APWR,
   /** Auto Increment Read Write */
   EC_CMD_APRW,
   /** Configured Address Read */
   EC_CMD_FPRD,
   /** Configured Address Write */
   EC_CMD_FPWR,
   /** Configured Address Read Write */
   EC_CMD_FPRW,
   /** Broadcast Read */
   EC_CMD_BRD,
   /** Broadcast Write */
   EC_CMD_BWR,
   /** Broadcast Read Write */
   EC_CMD_BRW,
   /** Logical Memory Read */
   EC_CMD_LRD,
   /** Logical Memory Write */
   EC_CMD_LWR,
   /** Logical Memory Read Write */
   EC_CMD_LRW,
   /** Auto Increment Read Multiple Write */
   EC_CMD_ARMW,
   /** Configured Read Multiple Write */
   EC_CMD_FRMW
   /** Reserved */
} ec_cmdtype;

/** Ethercat EEprom command types */
typedef enum
{
   /** No operation */
   EC_ECMD_NOP         = 0x0000,
   /** Read */
   EC_ECMD_READ        = 0x0100,
   /** Write */
   EC_ECMD_WRITE       = 0x0201,
   /** Reload */
   EC_ECMD_RELOAD      = 0x0300
} ec_ecmdtype;

/** EEprom state machine read size */
#define EC_ESTAT_R64    0x0040
/** EEprom state machine busy flag */
#define EC_ESTAT_BUSY   0x8000
/** EEprom state machine error flag mask */
#define EC_ESTAT_EMASK  0x7800
/** EEprom state machine error acknowledge */
#define EC_ESTAT_NACK   0x2000

/* Ethercat SSI (Slave Information Interface) */

/** Start address SII sections in Eeprom */
#define ECT_SII_START   0x0040

enum
{
   /** SII category strings */
   ECT_SII_STRING      = 10,
   /** SII category general */
   ECT_SII_GENERAL     = 30,
   /** SII category FMMU */
   ECT_SII_FMMU        = 40,
   /** SII category SM */
   ECT_SII_SM          = 41,
   /** SII category PDO */
   ECT_SII_PDO         = 50
};

/** Item offsets in SII general section */
enum
{
   ECT_SII_MANUF       = 0x0008,
   ECT_SII_ID          = 0x000a,
   ECT_SII_REV         = 0x000c,
   ECT_SII_BOOTRXMBX   = 0x0014,
   ECT_SII_BOOTTXMBX   = 0x0016,
   ECT_SII_MBXSIZE     = 0x0019,
   ECT_SII_TXMBXADR    = 0x001a,
   ECT_SII_RXMBXADR    = 0x0018,
   ECT_SII_MBXPROTO    = 0x001c
};

/** Mailbox types definitions */
enum
{
   /** Error mailbox type */
   ECT_MBXT_ERR        = 0x00,
   /** ADS over EtherCAT mailbox type */
   ECT_MBXT_AOE,
   /** Ethernet over EtherCAT mailbox type */
   ECT_MBXT_EOE,
   /** CANopen over EtherCAT mailbox type */
   ECT_MBXT_COE,
   /** File over EtherCAT mailbox type */
   ECT_MBXT_FOE,
   /** Servo over EtherCAT mailbox type */
   ECT_MBXT_SOE,
   /** Vendor over EtherCAT mailbox type */
   ECT_MBXT_VOE        = 0x0f
};

/** CoE mailbox types */
enum
{
   ECT_COES_EMERGENCY  = 0x01,
   ECT_COES_SDOREQ,
   ECT_COES_SDORES,
   ECT_COES_TXPDO,
   ECT_COES_RXPDO,
   ECT_COES_TXPDO_RR,
   ECT_COES_RXPDO_RR,
   ECT_COES_SDOINFO
};

/** CoE SDO commands */
enum
{
   ECT_SDO_DOWN_INIT    = 0x21,
   ECT_SDO_DOWN_EXP     = 0x23,
   ECT_SDO_DOWN_INIT_CA = 0x31,
   ECT_SDO_UP_REQ       = 0x40,
   ECT_SDO_UP_REQ_CA    = 0x50,
   ECT_SDO_SEG_UP_REQ   = 0x60,
   ECT_SDO_ABORT        = 0x80
};

/** CoE Object Description commands */
enum
{
   ECT_GET_ODLIST_REQ  = 0x01,
   ECT_GET_ODLIST_RES  = 0x02,
   ECT_GET_OD_REQ      = 0x03,
   ECT_GET_OD_RES      = 0x04,
   ECT_GET_OE_REQ      = 0x05,
   ECT_GET_OE_RES      = 0x06,
   ECT_SDOINFO_ERROR   = 0x07
};

/** FoE opcodes */
enum
{
   ECT_FOE_READ        = 0x01,
   ECT_FOE_WRITE,
   ECT_FOE_DATA,
   ECT_FOE_ACK,
   ECT_FOE_ERROR,
   ECT_FOE_BUSY
};

/** SoE opcodes */
enum
{
   ECT_SOE_READREQ     = 0x01,
   ECT_SOE_READRES,
   ECT_SOE_WRITEREQ,
   ECT_SOE_WRITERES,
   ECT_SOE_NOTIFICATION,
   ECT_SOE_EMERGENCY
};

/** Ethercat registers */
enum
{
   ECT_REG_TYPE        = 0x0000,
   ECT_REG_PORTDES     = 0x0007,
   ECT_REG_ESCSUP      = 0x0008,
   ECT_REG_STADR       = 0x0010,
   ECT_REG_ALIAS       = 0x0012,
   ECT_REG_DLCTL       = 0x0100,
   ECT_REG_DLPORT      = 0x0101,
   ECT_REG_DLALIAS     = 0x0103,
   ECT_REG_DLSTAT      = 0x0110,
   ECT_REG_ALCTL       = 0x0120,
   ECT_REG_ALSTAT      = 0x0130,
   ECT_REG_ALSTATCODE  = 0x0134,
   ECT_REG_PDICTL      = 0x0140,
   ECT_REG_IRQMASK     = 0x0200,
   ECT_REG_RXERR       = 0x0300,
   ECT_REG_FRXERR      = 0x0308,
   ECT_REG_EPUECNT     = 0x030C,
   ECT_REG_PECNT       = 0x030D,
   ECT_REG_PECODE      = 0x030E,
   ECT_REG_LLCNT       = 0x0310,
   ECT_REG_WDCNT       = 0x0442,
   ECT_REG_EEPCFG      = 0x0500,
   ECT_REG_EEPCTL      = 0x0502,
   ECT_REG_EEPSTAT     = 0x0502,
   ECT_REG_EEPADR      = 0x0504,
   ECT_REG_EEPDAT      = 0x0508,
   ECT_REG_FMMU0       = 0x0600,
   ECT_REG_FMMU1       = ECT_REG_FMMU0 + 0x10,
   ECT_REG_FMMU2       = ECT_REG_FMMU1 + 0x10,
   ECT_REG_FMMU3       = ECT_REG_FMMU2 + 0x10,
   ECT_REG_SM0         = 0x0800,
   ECT_REG_SM1         = ECT_REG_SM0 + 0x08,
   ECT_REG_SM2         = ECT_REG_SM1 + 0x08,
   ECT_REG_SM3         = ECT_REG_SM2 + 0x08,
   ECT_REG_SM0STAT     = ECT_REG_SM0 + 0x05,
   ECT_REG_SM1STAT     = ECT_REG_SM1 + 0x05,
   ECT_REG_SM1ACT      = ECT_REG_SM1 + 0x06,
   ECT_REG_SM1CONTR    = ECT_REG_SM1 + 0x07,
   ECT_REG_DCTIME0     = 0x0900,
   ECT_REG_DCTIME1     = 0x0904,
   ECT_REG_DCTIME2     = 0x0908,
   ECT_REG_DCTIME3     = 0x090C,
   ECT_REG_DCSYSTIME   = 0x0910,
   ECT_REG_DCSOF       = 0x0918,
   ECT_REG_DCSYSOFFSET = 0x0920,
   ECT_REG_DCSYSDELAY  = 0x0928,
   ECT_REG_DCSYSDIFF   = 0x092C,
   ECT_REG_DCSPEEDCNT  = 0x0930,
   ECT_REG_DCTIMEFILT  = 0x0934,
   ECT_REG_DCCUC       = 0x0980,
   ECT_REG_DCSYNCACT   = 0x0981,
   ECT_REG_DCSTART0    = 0x0990,
   ECT_REG_DCCYCLE0    = 0x09A0,
   ECT_REG_DCCYCLE1    = 0x09A4
};

/** standard SDO Sync Manager Communication Type */
#define ECT_SDO_SMCOMMTYPE      0x1c00
/** standard SDO PDO assignment */
#define ECT_SDO_PDOASSIGN       0x1c10
/** standard SDO RxPDO assignment */
#define ECT_SDO_RXPDOASSIGN     0x1c12
/** standard SDO TxPDO assignment */
#define ECT_SDO_TXPDOASSIGN     0x1c13

/** Ethercat packet type */
#define ETH_P_ECAT              0x88A4

/** Error types */
typedef enum
{
   EC_ERR_TYPE_SDO_ERROR         = 0,
   EC_ERR_TYPE_EMERGENCY         = 1,
   EC_ERR_TYPE_PACKET_ERROR      = 3,
   EC_ERR_TYPE_SDOINFO_ERROR     = 4,
   EC_ERR_TYPE_FOE_ERROR         = 5,
   EC_ERR_TYPE_FOE_BUF2SMALL     = 6,
   EC_ERR_TYPE_FOE_PACKETNUMBER  = 7,
   EC_ERR_TYPE_SOE_ERROR         = 8,
   EC_ERR_TYPE_MBX_ERROR         = 9,
   EC_ERR_TYPE_FOE_FILE_NOTFOUND = 10
} ec_err_type;

/** Struct to retrieve errors. */
typedef struct
{
   /** Time at which the error was generated. */
   ec_timet Time;
   /** Signal bit, error set but not read */
   boolean     Signal;
   /** Slave number that generated the error */
   uint16      Slave;
   /** CoE SDO index that generated the error */
   uint16      Index;
   /** CoE SDO subindex that generated the error */
   uint8       SubIdx;
   /** Type of error */
   ec_err_type Etype;
   union
   {
      /** General abortcode */
      int32   AbortCode;
      /** Specific error for Emergency mailbox */
      struct
      {
         uint16  ErrorCode;
         uint8   ErrorReg;
         uint8   b1;
         uint16  w1;
         uint16  w2;
      };
   };
} ec_errort;

/** Helper macros */
/** Macro to make a word from 2 bytes */
#define MK_WORD(msb, lsb)   ((((uint16)(msb))<<8) | (lsb))
/** Macro to get hi byte of a word */
#define HI_BYTE(w)          ((w) >> 8)
/** Macro to get low byte of a word */
#define LO_BYTE(w)          ((w) & 0x00ff)
/** Macro to swap hi and low byte of a word */
#define SWAP(w)             ((((w)& 0xff00) >> 8) | (((w) & 0x00ff) << 8))
/** Macro to get hi word of a dword */
#define LO_WORD(l)          ((l) & 0xffff)
/** Macro to get hi word of a dword */
#define HI_WORD(l)          ((l) >> 16)

#define get_unaligned(ptr) \
  ({ __typeof__(*(ptr)) __tmp; memcpy(&__tmp, (ptr), sizeof(*(ptr))); __tmp; })

#define put_unaligned32(val, ptr)        \
  (memcpy((ptr), &(val), 4))

#define put_unaligned64(val, ptr)        \
  (memcpy((ptr), &(val), 8))

#if !defined(EC_BIG_ENDIAN) && defined(EC_LITTLE_ENDIAN)

  #define htoes(A) (A)
  #define htoel(A) (A)
  #define htoell(A) (A)
  #define etohs(A) (A)
  #define etohl(A) (A)
  #define etohll(A) (A)

#elif !defined(EC_LITTLE_ENDIAN) && defined(EC_BIG_ENDIAN)

  #define htoes(A) ((((uint16)(A) & 0xff00) >> 8) | \
                    (((uint16)(A) & 0x00ff) << 8))
  #define htoel(A) ((((uint32)(A) & 0xff000000) >> 24) | \
                    (((uint32)(A) & 0x00ff0000) >> 8)  | \
                    (((uint32)(A) & 0x0000ff00) << 8)  | \
                    (((uint32)(A) & 0x000000ff) << 24))
  #define htoell(A) ((((uint64)(A) & (uint64)0xff00000000000000ULL) >> 56) | \
                     (((uint64)(A) & (uint64)0x00ff000000000000ULL) >> 40) | \
                     (((uint64)(A) & (uint64)0x0000ff0000000000ULL) >> 24) | \
                     (((uint64)(A) & (uint64)0x000000ff00000000ULL) >> 8)  | \
                     (((uint64)(A) & (uint64)0x00000000ff000000ULL) << 8)  | \
                     (((uint64)(A) & (uint64)0x0000000000ff0000ULL) << 24) | \
                     (((uint64)(A) & (uint64)0x000000000000ff00ULL) << 40) | \
                     (((uint64)(A) & (uint64)0x00000000000000ffULL) << 56))

  #define etohs  htoes
  #define etohl  htoel
  #define etohll htoell

#else

  #error "Must define one of EC_BIG_ENDIAN or EC_LITTLE_ENDIAN"

#endif

#ifdef __cplusplus
}
#endif

#endif /* _EC_TYPE_H */
