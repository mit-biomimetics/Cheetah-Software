#include <pthread.h>
#include <rt/rt_interface_lcm.h>

static lcm::LCM *g_lcm;

#include <rt/rt_sbus.h>

static pthread_mutex_t lcm_get_set_mutex =
    PTHREAD_MUTEX_INITIALIZER; /**< mutex to protect gui settings coming over
                                  LCM */

int iterations_since_last_lcm = 0;

/**
 * @brief      Increments the number of control iterations since the last LCM
 * packet
 */
void control_iteration_lcm() { iterations_since_last_lcm++; }

/**
 * @brief      Gets the iterations since the last lcm.
 *
 * @return     The iterations since last lcm.
 */
int get_iterations_since_last_lcm() { return iterations_since_last_lcm; }

// Controller Settings
#include <gui_main_control_settings_t.hpp>
// volatile gui_main_control_settings_t main_control_settings;
gui_main_control_settings_t main_control_settings;

#include <rc_channels_t.hpp>
// volatile rc_channels_t rc_channels;
rc_channels_t rc_channels;

/* ------------------------- HANDLERS ------------------------- */

// Controller Settings

void Handler::main_control_settings_handler(
    const lcm::ReceiveBuffer *rbuf, const std::string &channel,
    const gui_main_control_settings_t *msg) {
  (void)rbuf;
  (void)channel;

  iterations_since_last_lcm = 0;
  pthread_mutex_lock(&lcm_get_set_mutex);
  main_control_settings = *msg;
  pthread_mutex_unlock(&lcm_get_set_mutex);
}

void Handler::rc_channels_handler(const lcm::ReceiveBuffer *rbuf,
                                  const std::string &channel,
                                  const rc_channels_t *msg) {
  (void)rbuf;
  (void)channel;

  pthread_mutex_lock(&lcm_get_set_mutex);
  rc_channels = *msg;
  pthread_mutex_unlock(&lcm_get_set_mutex);
}
void get_main_control_settings(void *settings) {
  pthread_mutex_lock(&lcm_get_set_mutex);
  v_memcpy(settings, &main_control_settings, sizeof(main_control_settings));
  pthread_mutex_unlock(&lcm_get_set_mutex);
}

void get_rc_channels(void *settings) {
  pthread_mutex_lock(&lcm_get_set_mutex);
  v_memcpy(settings, &rc_channels, sizeof(rc_channels));
  pthread_mutex_unlock(&lcm_get_set_mutex);
}

/**
 * @brief      Function which handles the completion of an SBUS Packet and
 * overrides the LCM control settings as desired.
 */
void sbus_packet_complete() {
  pthread_mutex_lock(&lcm_get_set_mutex);

  int ch1 = read_sbus_channel(0);
  int ch2 = read_sbus_channel(1);
  int ch3 = read_sbus_channel(2);
  int ch4 = read_sbus_channel(3);
  int ch7 = read_sbus_channel(6);
  int ch9 = read_sbus_channel(8);
  int ch10 = read_sbus_channel(9);
  int ch11 = read_sbus_channel(10);
  int ch13 = read_sbus_channel(12);

  // velocity scales between 1.0 and 3.0
  float v_scale = 2.0f * (((float)(ch7 - 172)) / 1811.0f) + 1.0f;
  // Ignore commands if switched
  if (ch11 != 1811) {
    // oh shit switch
    if (ch10 == 172) {
      main_control_settings.mode = 0;
    } else if (ch10 == 992) {
      main_control_settings.mode = 12;
    } else if (ch10 == 1811) {
      main_control_settings.mode = 11;
    }

    // Use the joysticks for velocity and yaw control in locomotion gaits
    if (main_control_settings.mode == 11) {
      // Analog channels return a value between ~200 and 1800, with sticks
      // centered at 1000

      if (ch9 == 1811) {
        main_control_settings.variable[0] = 4;
      }  // Stand
      else if (ch9 == 992) {
        main_control_settings.variable[0] = 0;
      }  // Trot
      else if ((ch9 == 172) && (ch13 == 1811)) {
        main_control_settings.variable[0] = 5;
      }  // Trot run
      else if ((ch9 == 172) && (ch13 == 992)) {
        main_control_settings.variable[0] = 2;
      }  // Pronk
      else if ((ch9 == 172) && (ch13 == 172)) {
        main_control_settings.variable[0] = 1;
      }  // Bound

      // If using the stand gait, control orientations
      if (main_control_settings.variable[0] == 4) {
        main_control_settings.rpy_des[0] = ((float)ch4 - 1000) * .0003f;
        main_control_settings.rpy_des[1] = ((float)ch3 - 1000) * .0005f;
        main_control_settings.rpy_des[2] = ((float)ch2 - 1000) * .001f;

        main_control_settings.p_des[0] = 0;
        main_control_settings.p_des[1] = 0;
        main_control_settings.p_des[2] = 0.21 + ((float)ch1 - 1000) * .0001f;
      }
      // For all other gaits, control the velocities
      else {
        main_control_settings.v_des[0] = v_scale * ((float)ch1 - 1000) * .002f;
        main_control_settings.v_des[1] =
            -v_scale * ((float)ch4 - 1000) * .0005f;
        main_control_settings.v_des[2] = 0;
        main_control_settings.p_des[2] = 0.21;

        main_control_settings.omega_des[0] = 0;
        main_control_settings.omega_des[1] = 0;
        main_control_settings.omega_des[2] =
            -v_scale * ((float)ch2 - 1000) * .003f;
      }
    }
    // For standing modes (mpc or balance qp) control orientations
    else if (main_control_settings.mode == 10 ||
             main_control_settings.mode == 3) {
      main_control_settings.rpy_des[0] = ((float)ch4 - 1000) * .0003f;
      main_control_settings.rpy_des[1] = ((float)ch3 - 1000) * .0005f;
      main_control_settings.rpy_des[2] = ((float)ch2 - 1000) * .001f;

      main_control_settings.p_des[0] = 0;
      main_control_settings.p_des[1] = 0;
      main_control_settings.p_des[2] = 0.21 + ((float)ch1 - 1000) * .0001f;
    }
  }
  // Use the joysticks for orientation and height control in standing mode

  pthread_mutex_unlock(&lcm_get_set_mutex);
  // printf("[RT Interface LCM] Got SBUS Packet\n");
}

/**
 * @brief      Initializer for all handlers related to the interface LCM streams
 *
 * @param      main_lcm  A pointer to the main lcm object used in the lcm thread
 */
void init_interface_lcm(lcm::LCM *main_lcm) {
  printf("[RT Interface LCM] Initializing...\n");
  g_lcm = main_lcm;

  Handler handlerObj;

  g_lcm->subscribe("INTERFACE_gui_main_control_settings",
                   &Handler::main_control_settings_handler, &handlerObj);
  g_lcm->subscribe("INTERFACE_rc_channels", &Handler::rc_channels_handler,
                   &handlerObj);

  main_control_settings.enable = 1;

  printf("[RT Interface LCM] Done\n");

  // printf("Initial PFOOT SETTINGS %lf\n",
  // state_estimator_settings.process_noise_pfoot);
}

void *v_memcpy(void *dest, volatile void *src, size_t n) {
  void *src_2 = (void *)src;
  return memcpy(dest, src_2, n);
}
