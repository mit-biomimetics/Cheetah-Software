/*
 * Simple Open EtherCAT Master Library
 *
 * File    : oshw.c
 * Version : 1.3.1
 * Date    : 11-03-2015
 * Copyright (C) 2005-2015 Speciaal Machinefabriek Ketels v.o.f.
 * Copyright (C) 2005-2015 Arthur Ketels
 * Copyright (C) 2008-2009 TU/e Technische Universiteit Eindhoven
 * Copyright (C) 2012-2015 rt-labs AB , Sweden
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

#include <sys/ioctl.h>
#include <net/if.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "oshw.h"

/**
 * Host to Network byte order (i.e. to big endian).
 *
 * Note that Ethercat uses little endian byte order, except for the Ethernet
 * header which is big endian as usual.
 */
uint16 oshw_htons(uint16 host)
{
   uint16 network = htons (host);
   return network;
}

/**
 * Network (i.e. big endian) to Host byte order.
 *
 * Note that Ethercat uses little endian byte order, except for the Ethernet
 * header which is big endian as usual.
 */
uint16 oshw_ntohs(uint16 network)
{
   uint16 host = ntohs (network);
   return host;
}

/** Create list over available network adapters.
 * @return First element in linked list of adapters
 */
ec_adaptert * oshw_find_adapters(void)
{
   int i;
   int string_len;
   struct if_nameindex *ids;
   ec_adaptert * adapter;
   ec_adaptert * prev_adapter;
   ec_adaptert * ret_adapter = NULL;


   /* Iterate all devices and create a local copy holding the name and
    * decsription.
    */

   ids = if_nameindex ();
   for(i = 0; ids[i].if_index != 0; i++)
   {
      adapter = (ec_adaptert *)malloc(sizeof(ec_adaptert));
      /* If we got more than one adapter save link list pointer to previous
       * adapter.
       * Else save as pointer to return.
       */
      if (i)
      {
         prev_adapter->next = adapter;
      }
      else
      {
         ret_adapter = adapter;
      }

      /* fetch description and name, in linux we use the same on both */
      adapter->next = NULL;

      if (ids[i].if_name)
      {
          string_len = strlen(ids[i].if_name);
          if (string_len > (EC_MAXLEN_ADAPTERNAME - 1))
          {
             string_len = EC_MAXLEN_ADAPTERNAME - 1;
          }
          strncpy(adapter->name, ids[i].if_name,string_len);
          adapter->name[string_len] = '\0';
          strncpy(adapter->desc, ids[i].if_name,string_len);
          adapter->desc[string_len] = '\0';
      }
      else
      {
         adapter->name[0] = '\0';
         adapter->desc[0] = '\0';
      }

      prev_adapter = adapter;
   }

   if_freenameindex (ids);

   return ret_adapter;
}

/** Free memory allocated memory used by adapter collection.
 * @param[in] adapter = First element in linked list of adapters
 * EC_NOFRAME.
 */
void oshw_free_adapters(ec_adaptert * adapter)
{
   ec_adaptert * next_adapter;
   /* Iterate the linked list and free all elemnts holding
    * adapter information
    */
   if(adapter)
   {
      next_adapter = adapter->next;
      free (adapter);
      while (next_adapter)
      {
         adapter = next_adapter;
         next_adapter = adapter->next;
         free (adapter);
      }
   }
}
