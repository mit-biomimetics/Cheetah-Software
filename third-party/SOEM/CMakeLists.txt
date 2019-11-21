include_directories(osal)
include_directories(osal/linux)
include_directories(oshw)
include_directories(oshw/linux)
include_directories(soem)

add_library(soem SHARED
        osal/linux/osal.c
        oshw/linux/nicdrv.c
        oshw/linux/oshw.c
        soem/ethercatbase.c
        soem/ethercatcoe.c
        soem/ethercatconfig.c
        soem/ethercatdc.c
        soem/ethercatfoe.c
        soem/ethercatmain.c
        soem/ethercatprint.c
        soem/ethercatsoe.c
        )