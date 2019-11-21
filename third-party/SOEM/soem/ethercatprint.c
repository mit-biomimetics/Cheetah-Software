/*
 * Simple Open EtherCAT Master Library
 *
 * File    : ethercatprint.c
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
 * Module to convert EtherCAT errors to readable messages.
 *
 * SDO abort messages and AL status codes are used to relay slave errors to
 * the user application. This module converts the binary codes to readble text.
 * For the defined error codes see the EtherCAT protocol documentation.
 */

#include <stdio.h>
#include "oshw.h"
#include "ethercattype.h"
#include "ethercatmain.h"

#define EC_MAXERRORNAME 127

/** SDO error list type definition */
typedef struct
{
   /** Error code returned from SDO */
   uint32        errorcode;
   /** Readable error description */
   char          errordescription[EC_MAXERRORNAME + 1];
} ec_sdoerrorlist_t;

/** AL status code list type definition */
typedef struct
{
   /** AL status code */
   uint16        ALstatuscode;
   /** Readable description */
   char          ALstatuscodedescription[EC_MAXERRORNAME + 1];
} ec_ALstatuscodelist_t;

/** SoE error list type definition */
typedef struct
{
   /** SoE error code */
   uint16        errorcode;
   /** Readable description */
   char          errordescription[EC_MAXERRORNAME + 1];
} ec_soeerrorlist_t;

/** MBX error list type definition */
typedef struct
{
   /** MBX error code */
   uint16              errorcode;
   /** Readable description */
   char                errordescription[EC_MAXERRORNAME + 1];
} ec_mbxerrorlist_t;

char estring[EC_MAXERRORNAME];

/** SDO error list definition */
const ec_sdoerrorlist_t ec_sdoerrorlist[] = {
   {0x00000000, "No error" },
   {0x05030000, "Toggle bit not changed" },
   {0x05040000, "SDO protocol timeout" },
   {0x05040001, "Client/Server command specifier not valid or unknown" },
   {0x05040005, "Out of memory" },
   {0x06010000, "Unsupported access to an object" },
   {0x06010001, "Attempt to read to a write only object" },
   {0x06010002, "Attempt to write to a read only object" },
   {0x06010003, "Subindex can not be written, SI0 must be 0 for write access" },
   {0x06010004, "SDO Complete access not supported for variable length objects" },
   {0x06010005, "Object length exceeds mailbox size" },
   {0x06010006, "Object mapped to RxPDO, SDO download blocked" },
   {0x06020000, "The object does not exist in the object directory" },
   {0x06040041, "The object can not be mapped into the PDO" },
   {0x06040042, "The number and length of the objects to be mapped would exceed the PDO length" },
   {0x06040043, "General parameter incompatibility reason" },
   {0x06040047, "General internal incompatibility in the device" },
   {0x06060000, "Access failed due to a hardware error" },
   {0x06070010, "Data type does not match, length of service parameter does not match" },
   {0x06070012, "Data type does not match, length of service parameter too high" },
   {0x06070013, "Data type does not match, length of service parameter too low" },
   {0x06090011, "Subindex does not exist" },
   {0x06090030, "Value range of parameter exceeded (only for write access)" },
   {0x06090031, "Value of parameter written too high" },
   {0x06090032, "Value of parameter written too low" },
   {0x06090036, "Maximum value is less than minimum value" },
   {0x08000000, "General error" },
   {0x08000020, "Data cannot be transferred or stored to the application" },
   {0x08000021, "Data cannot be transferred or stored to the application because of local control" },
   {0x08000022, "Data cannot be transferred or stored to the application because of the present device state" },
   {0x08000023, "Object dictionary dynamic generation fails or no object dictionary is present" },
   {0xffffffff, "Unknown" }
};

/** AL status code list definition */
const ec_ALstatuscodelist_t ec_ALstatuscodelist[] = {
   {0x0000 , "No error" },
   {0x0001 , "Unspecified error" },
   {0x0002 , "No memory" },
   {0x0011 , "Invalid requested state change" },
   {0x0012 , "Unknown requested state" },
   {0x0013 , "Bootstrap not supported" },
   {0x0014 , "No valid firmware" },
   {0x0015 , "Invalid mailbox configuration" },
   {0x0016 , "Invalid mailbox configuration" },
   {0x0017 , "Invalid sync manager configuration" },
   {0x0018 , "No valid inputs available" },
   {0x0019 , "No valid outputs" },
   {0x001A , "Synchronization error" },
   {0x001B , "Sync manager watchdog" },
   {0x001C , "Invalid sync Manager types" },
   {0x001D , "Invalid output configuration" },
   {0x001E , "Invalid input configuration" },
   {0x001F , "Invalid watchdog configuration" },
   {0x0020 , "Slave needs cold start" },
   {0x0021 , "Slave needs INIT" },
   {0x0022 , "Slave needs PREOP" },
   {0x0023 , "Slave needs SAFEOP" },
   {0x0024 , "Invalid input mapping" },
   {0x0025 , "Invalid output mapping" },
   {0x0026 , "Inconsistent settings" },
   {0x0027 , "Freerun not supported" },
   {0x0028 , "Synchronisation not supported" },
   {0x0029 , "Freerun needs 3buffer mode" },
   {0x002A , "Background watchdog" },
   {0x002B , "No valid Inputs and Outputs" },
   {0x002C , "Fatal sync error" },
   {0x002D , "No sync error" }, // was "Invalid Output FMMU Configuration"
   {0x002E , "Invalid input FMMU configuration" },
   {0x0030 , "Invalid DC SYNC configuration" },
   {0x0031 , "Invalid DC latch configuration" },
   {0x0032 , "PLL error" },
   {0x0033 , "DC sync IO error" },
   {0x0034 , "DC sync timeout error" },
   {0x0035 , "DC invalid sync cycle time" },
   {0x0036 , "DC invalid sync0 cycle time" },
   {0x0037 , "DC invalid sync1 cycle time" },
   {0x0041 , "MBX_AOE" },
   {0x0042 , "MBX_EOE" },
   {0x0043 , "MBX_COE" },
   {0x0044 , "MBX_FOE" },
   {0x0045 , "MBX_SOE" },
   {0x004F , "MBX_VOE" },
   {0x0050 , "EEPROM no access" },
   {0x0051 , "EEPROM error" },
   {0x0060 , "Slave restarted locally" },
   {0x0061 , "Device identification value updated" },
   {0x00f0 , "Application controller available" },
   {0xffff , "Unknown" }
};

/** SoE error list definition */
const ec_soeerrorlist_t ec_soeerrorlist[] = {
   {0x0000, "No error" },
   {0x1001, "No IDN" },
   {0x1009, "Invalid access to element 1" },
   {0x2001, "No Name" },
   {0x2002, "Name transmission too short" },
   {0x2003, "Name transmission too long" },
   {0x2004, "Name cannot be changed (read only)" },
   {0x2005, "Name is write-protected at this time" },
   {0x3002, "Attribute transmission too short" },
   {0x3003, "Attribute transmission too long" },
   {0x3004, "Attribute cannot be changed (read only)" },
   {0x3005, "Attribute is write-protected at this time" },
   {0x4001, "No units" },
   {0x4002, "Unit transmission too short" },
   {0x4003, "Unit transmission too long" },
   {0x4004, "Unit cannot be changed (read only)" },
   {0x4005, "Unit is write-protected at this time" },
   {0x5001, "No minimum input value" },
   {0x5002, "Minimum input value transmission too short" },
   {0x5003, "Minimum input value transmission too long" },
   {0x5004, "Minimum input value cannot be changed (read only)" },
   {0x5005, "Minimum input value is write-protected at this time" },
   {0x6001, "No maximum input value" },
   {0x6002, "Maximum input value transmission too short" },
   {0x6003, "Maximum input value transmission too long" },
   {0x6004, "Maximum input value cannot be changed (read only)" },
   {0x6005, "Maximum input value is write-protected at this time" },
   {0x7002, "Operation data transmission too short" },
   {0x7003, "Operation data transmission too long" },
   {0x7004, "Operation data cannot be changed (read only)" },
   {0x7005, "Operation data is write-protected at this time (state)" },
   {0x7006, "Operation data is smaller than the minimum input value" },
   {0x7007, "Operation data is smaller than the maximum input value" },
   {0x7008, "Invalid operation data:Configured IDN will not be supported" },
   {0x7009, "Operation data write protected by a password" },
   {0x700A, "Operation data is write protected, it is configured cyclically" },
   {0x700B, "Invalid indirect addressing: (e.g., data container, list handling)" },
   {0x700C, "Operation data is write protected, due to other settings" },
   {0x700D, "Reserved" },
   {0x7010, "Procedure command already active" },
   {0x7011, "Procedure command not interruptible" },
   {0x7012, "Procedure command at this time not executable (state)" },
   {0x7013, "Procedure command not executable (invalid or false parameters)" },
   {0x7014, "No data state" },
   {0x8001, "No default value" },
   {0x8002, "Default value transmission too long" },
   {0x8004, "Default value cannot be changed, read only" },
   {0x800A, "Invalid drive number" },
   {0x800A, "General error" },
   {0x800A, "No element addressed" },
   {0xffff, "Unknown" }
};

/** MBX error list definition */
const ec_mbxerrorlist_t ec_mbxerrorlist[] = {
   {0x0000, "No error" },
   {0x0001, "Syntax of 6 octet Mailbox Header is wrong" },
   {0x0002, "The mailbox protocol is not supported" },
   {0x0003, "Channel Field contains wrong value"},
   {0x0004, "The service is no supported"},
   {0x0005, "Invalid mailbox header"},
   {0x0006, "Length of received mailbox data is too short"},
   {0x0007, "No more memory in slave"},
   {0x0008, "The lenght of data is inconsistent"},
   {0xffff, "Unknown"}
};

/** Look up text string that belongs to SDO error code.
 *
 * @param[in] sdoerrorcode   = SDO error code as defined in EtherCAT protocol
 * @return readable string
 */
const char* ec_sdoerror2string( uint32 sdoerrorcode)
{
   int i = 0;

   while ( (ec_sdoerrorlist[i].errorcode != 0xffffffffUL) &&
           (ec_sdoerrorlist[i].errorcode != sdoerrorcode) )
   {
      i++;
   }

   return ec_sdoerrorlist[i].errordescription;
}

/** Look up text string that belongs to AL status code.
 *
 * @param[in] ALstatuscode   = AL status code as defined in EtherCAT protocol
 * @return readable string
 */
char* ec_ALstatuscode2string( uint16 ALstatuscode)
{
   int i = 0;

   while ( (ec_ALstatuscodelist[i].ALstatuscode != 0xffff) &&
           (ec_ALstatuscodelist[i].ALstatuscode != ALstatuscode) )
   {
      i++;
   }

   return (char *) ec_ALstatuscodelist[i].ALstatuscodedescription;
}

/** Look up text string that belongs to SoE error code.
 *
 * @param[in] errorcode   = SoE error code as defined in EtherCAT protocol
 * @return readable string
 */
char* ec_soeerror2string( uint16 errorcode)
{
   int i = 0;

   while ( (ec_soeerrorlist[i].errorcode != 0xffff) &&
           (ec_soeerrorlist[i].errorcode != errorcode) )
   {
      i++;
   }

   return (char *) ec_soeerrorlist[i].errordescription;
}

/** Look up text string that belongs to MBX error code.
 *
 * @param[in] errorcode   = MBX error code as defined in EtherCAT protocol
 * @return readable string
 */
char* ec_mbxerror2string( uint16 errorcode)
{
   int i = 0;

   while ( (ec_mbxerrorlist[i].errorcode != 0xffff) &&
           (ec_mbxerrorlist[i].errorcode != errorcode) )
   {
      i++;
   }

   return (char *) ec_mbxerrorlist[i].errordescription;
}

/** Look up error in ec_errorlist and convert to text string.
 *
 * @param[in]  context        = context struct
 * @return readable string
 */
char* ecx_elist2string(ecx_contextt *context)
{
   ec_errort Ec;
   char timestr[20];

   if (ecx_poperror(context, &Ec))
   {
      sprintf(timestr, "Time:%12.3f", Ec.Time.sec + (Ec.Time.usec / 1000000.0) );
      switch (Ec.Etype)
      {
         case EC_ERR_TYPE_SDO_ERROR:
         {
            sprintf(estring, "%s SDO slave:%d index:%4.4x.%2.2x error:%8.8x %s\n",
                    timestr, Ec.Slave, Ec.Index, Ec.SubIdx, (unsigned)Ec.AbortCode, ec_sdoerror2string(Ec.AbortCode));
            break;
         }
         case EC_ERR_TYPE_EMERGENCY:
         {
            sprintf(estring, "%s EMERGENCY slave:%d error:%4.4x\n",
                    timestr, Ec.Slave, Ec.ErrorCode);
            break;
         }
         case EC_ERR_TYPE_PACKET_ERROR:
         {
            sprintf(estring, "%s PACKET slave:%d index:%4.4x.%2.2x error:%d\n",
                    timestr, Ec.Slave, Ec.Index, Ec.SubIdx, Ec.ErrorCode);
            break;
         }
         case EC_ERR_TYPE_SDOINFO_ERROR:
         {
            sprintf(estring, "%s SDO slave:%d index:%4.4x.%2.2x error:%8.8x %s\n",
                    timestr, Ec.Slave, Ec.Index, Ec.SubIdx, (unsigned)Ec.AbortCode, ec_sdoerror2string(Ec.AbortCode));
            break;
         }
         case EC_ERR_TYPE_SOE_ERROR:
         {
            sprintf(estring, "%s SoE slave:%d IDN:%4.4x error:%4.4x %s\n",
                    timestr, Ec.Slave, Ec.Index, (unsigned)Ec.AbortCode, ec_soeerror2string(Ec.ErrorCode));
            break;
         }
         case EC_ERR_TYPE_MBX_ERROR:
         {
            sprintf(estring, "%s MBX slave:%d error:%4.4x %s\n",
                    timestr, Ec.Slave, Ec.ErrorCode, ec_mbxerror2string(Ec.ErrorCode));
            break;
         }
         default:
         {
            sprintf(estring, "%s error:%8.8x\n",
                    timestr, (unsigned)Ec.AbortCode);
            break;
         }
      }
      return (char*) estring;
   }
   else
   {
      return "";
   }
}

#ifdef EC_VER1
char* ec_elist2string(void)
{
   return ecx_elist2string(&ecx_context);
}
#endif
