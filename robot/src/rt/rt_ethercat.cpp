#include "rt/rt_ethercat.h"
#include <stdio.h>
#include <string.h>
#include <inttypes.h> // what is this?
#include <stdlib.h>
#include <time.h>
#include <mutex>
#include <SOEM/soem/ethercat.h>
#include "SOEM/soem/ethercat.h"
#include "SOEM/osal/osal.h"
#include "SOEM/osal/linux/osal_defs.h"
#include "rt/rt_ethercat.h"
#include "SimUtilities/ti_boardcontrol.h"


#define EC_TIMEOUTMON 500

static char IOmap[4096];
static OSAL_THREAD_HANDLE thread1;
static int expectedWKC;
static boolean needlf;
static volatile int wkc;
static boolean inOP;
static uint8 currentgroup = 0;

//#define ADAPTER_NAME "eno1"
#define ADAPTER_NAME "enp2s0"

static void degraded_handler() {
  //shut of gpio enables
  // estop();
  printf("[EtherCAT Error] Logging error...\n");
  time_t current_time = time(NULL);
  char* time_str = ctime(&current_time);
  printf("ESTOP. EtherCAT became degraded at %s.\n", time_str);
  printf("[EtherCAT Error] Stopping RT process.\n");
  exit(0);
}

static int run_ethercat(const char *ifname) {
  int i, oloop, iloop, chk;
  needlf = FALSE;
  inOP = FALSE;


  /* initialise SOEM, bind socket to ifname */
  if (ec_init(ifname))
  {
    printf("[EtherCAT Init] Initialization on device %s succeeded.\n",ifname);
    /* find and auto-config slaves */

    if ( ec_config_init(FALSE) > 0 )
    {
      printf("[EtherCAT Init] %d slaves found and configured.\n",ec_slavecount);
      if(ec_slavecount < 4)
      {
        printf("[RT EtherCAT] Warning: Expected %d legs, found %d.\n", 4, ec_slavecount);
      }

      ec_config_map(&IOmap);
      ec_configdc();

      printf("[EtherCAT Init] Mapped slaves.\n");
      /* wait for all slaves to reach SAFE_OP state */
      ec_statecheck(0, EC_STATE_SAFE_OP,  EC_TIMEOUTSTATE * 4);

      for(int slave_idx = 0; slave_idx < ec_slavecount; slave_idx++) {
        printf("[SLAVE %d]\n", slave_idx);
        printf("  IN  %d bytes, %d bits\n", ec_slave[slave_idx].Ibytes, ec_slave[slave_idx].Ibits);
        printf("  OUT %d bytes, %d bits\n", ec_slave[slave_idx].Obytes, ec_slave[slave_idx].Obits);
        printf("\n");
     }

      oloop = ec_slave[0].Obytes;
      if ((oloop == 0) && (ec_slave[0].Obits > 0)) oloop = 1;
      if (oloop > 8) oloop = 8;
      iloop = ec_slave[0].Ibytes;
      if ((iloop == 0) && (ec_slave[0].Ibits > 0)) iloop = 1;
      if (iloop > 8) iloop = 8;

      printf("[EtherCAT Init] segments : %d : %d %d %d %d\n",ec_group[0].nsegments ,ec_group[0].IOsegment[0],ec_group[0].IOsegment[1],ec_group[0].IOsegment[2],ec_group[0].IOsegment[3]);

      printf("[EtherCAT Init] Requesting operational state for all slaves...\n");
      expectedWKC = (ec_group[0].outputsWKC * 2) + ec_group[0].inputsWKC;
      printf("[EtherCAT Init] Calculated workcounter %d\n", expectedWKC);
      ec_slave[0].state = EC_STATE_OPERATIONAL;
      /* send one valid process data to make outputs in slaves happy*/
      ec_send_processdata();
      ec_receive_processdata(EC_TIMEOUTRET);
      /* request OP state for all slaves */
      ec_writestate(0);
      chk = 40;
      /* wait for all slaves to reach OP state */
      do
      {
        ec_send_processdata();
        ec_receive_processdata(EC_TIMEOUTRET);
        ec_statecheck(0, EC_STATE_OPERATIONAL, 50000);
      }
      while (chk-- && (ec_slave[0].state != EC_STATE_OPERATIONAL));

      if (ec_slave[0].state == EC_STATE_OPERATIONAL )
      {
        printf("[EtherCAT Init] Operational state reached for all slaves.\n");
        inOP = TRUE;
        return 1;

      }
      else
      {
        printf("[EtherCAT Error] Not all slaves reached operational state.\n");
        ec_readstate();

        for(i = 1; i<=ec_slavecount ; i++)
        {
          if(ec_slave[i].state != EC_STATE_OPERATIONAL)
          {

            printf("[EtherCAT Error] Slave %d State=0x%2.2x StatusCode=0x%4.4x : %s\n",
                    i, ec_slave[i].state, ec_slave[i].ALstatuscode, ec_ALstatuscode2string(ec_slave[i].ALstatuscode));
          }
        }
      }
    }
    else
    {
      printf("[EtherCAT Error] No slaves found!\n");
    }
  }
  else
  {
    printf("[EtherCAT Error] No socket connection on %s - are you running run.sh?\n",ifname);
  }
  return 0;
}

static int err_count = 0;
static int err_iteration_count = 0;
/**@brief EtherCAT errors are measured over this period of loop iterations */
#define K_ETHERCAT_ERR_PERIOD 100

/**@brief Maximum number of etherCAT errors before a fault per period of loop iterations */
#define K_ETHERCAT_ERR_MAX 20

static OSAL_THREAD_FUNC ecatcheck( void *ptr )
{
  (void)ptr;
  int slave = 0;
  while(1)
  {
    //count errors
    if(err_iteration_count > K_ETHERCAT_ERR_PERIOD)
    {
      err_iteration_count = 0;
      err_count = 0;
    }

    if(err_count > K_ETHERCAT_ERR_MAX)
    {
      //possibly shut down
      printf("[EtherCAT Error] EtherCAT connection degraded.\n");
      printf("[Simulink-Linux] Shutting down....\n");
      degraded_handler();
      break;
    }
    err_iteration_count++;

    if( inOP && ((wkc < expectedWKC) || ec_group[currentgroup].docheckstate))
    {
      if (needlf)
      {
        needlf = FALSE;
        printf("\n");
      }
      /* one ore more slaves are not responding */
      ec_group[currentgroup].docheckstate = FALSE;
      ec_readstate();
      for (slave = 1; slave <= ec_slavecount; slave++)
      {
        if ((ec_slave[slave].group == currentgroup) && (ec_slave[slave].state != EC_STATE_OPERATIONAL))
        {
          ec_group[currentgroup].docheckstate = TRUE;
          if (ec_slave[slave].state == (EC_STATE_SAFE_OP + EC_STATE_ERROR))
          {
            printf("[EtherCAT Error] Slave %d is in SAFE_OP + ERROR, attempting ack.\n", slave);
            ec_slave[slave].state = (EC_STATE_SAFE_OP + EC_STATE_ACK);
            ec_writestate(slave);
            err_count++;
          }
          else if(ec_slave[slave].state == EC_STATE_SAFE_OP)
          {
            printf("[EtherCAT Error] Slave %d is in SAFE_OP, change to OPERATIONAL.\n", slave);
            ec_slave[slave].state = EC_STATE_OPERATIONAL;
            ec_writestate(slave);
            err_count++;
          }
          else if(ec_slave[slave].state > 0)
          {
            if (ec_reconfig_slave(slave, EC_TIMEOUTMON))
            {
              ec_slave[slave].islost = FALSE;
              printf("[EtherCAT Status] Slave %d reconfigured\n",slave);
            }
          }
          else if(!ec_slave[slave].islost)
          {
            /* re-check state */
            ec_statecheck(slave, EC_STATE_OPERATIONAL, EC_TIMEOUTRET);
            if (!ec_slave[slave].state)
            {
              ec_slave[slave].islost = TRUE;
              printf("[EtherCAT Error] Slave %d lost\n",slave);
              err_count++;
            }
          }
        }
        if (ec_slave[slave].islost)
        {
          if(!ec_slave[slave].state)
          {
            if (ec_recover_slave(slave, EC_TIMEOUTMON))
            {
              ec_slave[slave].islost = FALSE;
              printf("[EtherCAT Status] Slave %d recovered\n",slave);
            }
          }
          else
          {
            ec_slave[slave].islost = FALSE;
            printf("[EtherCAT Status] Slave %d found\n",slave);
          }
        }
      }
      if(!ec_group[currentgroup].docheckstate)
        printf("[EtherCAT Status] All slaves resumed OPERATIONAL.\n");
    }
    osal_usleep(50000);
  }
}

void rt_ethercat_init()
{

  printf("[EtherCAT] Initializing EtherCAT\n");
  //initialize monitoring thread
  osal_thread_create((void*)&thread1, 128000, (void*)&ecatcheck, (void*) &ctime);

  //try initialization until it succeeds
  int i;
  int rc;
  for(i = 1; i < 100; i++)
  {
    printf("[EtherCAT] Attempting to start EtherCAT, try %d of 100.\n", i);
    rc = run_ethercat(ADAPTER_NAME); // todo?
    if(rc) break;
    osal_usleep(1000000);
  }
  if(rc) printf("[EtherCAT] EtherCAT successfully initialized on attempt %d \n", i);
  else
  {
    printf("[EtherCAT Error] Failed to initialize EtherCAT after 100 tries. \n");
  }
  
}

static int wkc_err_count = 0;
static int wkc_err_iteration_count = 0;

static std::mutex command_mutex, data_mutex;

//initiate etherCAT communication
/** @brief Send and receive data over EtherCAT
 *
 * In Simulation, send data over LCM
 * On the robt, verify the EtherCAT connection is still healthy, send data, receive data, and check for lost packets
 */
void rt_ethercat_run()
{
  //check connection
  if(wkc_err_iteration_count > K_ETHERCAT_ERR_PERIOD)
  {
    wkc_err_count = 0;
    wkc_err_iteration_count = 0;
  }
  if(wkc_err_count > K_ETHERCAT_ERR_MAX)
  {
    printf("[EtherCAT Error] Error count too high!\n");
    //program terminates in degraded handler.
    degraded_handler();
  }

  //send
  command_mutex.lock();
  ec_send_processdata();
  command_mutex.unlock();

  //receive
  data_mutex.lock();
  wkc = ec_receive_processdata(EC_TIMEOUTRET);
  data_mutex.unlock();

  //check for dropped packet
  if(wkc < expectedWKC)
  {
    printf("\x1b[31m[EtherCAT Error] Dropped packet (Bad WKC!)\x1b[0m\n");
    wkc_err_count++;
  }
  wkc_err_iteration_count++;
}

void rt_ethercat_get_data(TiBoardData* data) {
  data_mutex.lock();
  for(int slave = 0; slave < 4; slave++) {
    TiBoardData* slave_src = (TiBoardData*)(ec_slave[slave + 1].inputs);
    if(slave_src)
      data[slave] = *(TiBoardData*)(ec_slave[slave + 1].inputs);
  }
  data_mutex.unlock();

}
void rt_ethercat_set_command(TiBoardCommand* command) {
  command_mutex.lock();
  for(int slave = 0; slave < 4; slave++) {
    TiBoardCommand* slave_dest = (TiBoardCommand*)(ec_slave[slave + 1].outputs);
    if(slave_dest)
      *(TiBoardCommand*)(ec_slave[slave + 1].outputs) = command[slave];
  }
  command_mutex.unlock();
}
