#include <pthread.h>
#include <rt/rt_interface_lcm.h>
#include "Utilities/EdgeTrigger.h"

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

EdgeTrigger<int> mode_edge_trigger(0);
EdgeTrigger<TaranisSwitchState> backflip_prep_edge_trigger(SWITCH_UP);
TaranisSwitchState initial_backflip_go_switch = SWITCH_DOWN;

void sbus_packet_complete() {
  Taranis_X7_data data;
  update_taranis_x7(&data);

  float v_scale = data.knobs[0]*1.2f + 1.2f; // from 0.5 to 1.5

  auto estop_switch = data.right_lower_right_switch;
  auto backflip_prepare_switch = data.left_lower_left_switch;
  auto backflip_go_switch = data.left_upper_switch;

  auto left_select = data.left_lower_right_switch;
  auto right_select = data.right_lower_left_switch;

  int selected_mode = 0;

  switch(estop_switch) {

    case SWITCH_UP: // ESTOP
      selected_mode = RC_mode::OFF;
      break;

    case SWITCH_MIDDLE: // recover
      selected_mode = RC_mode::RECOVERY_STAND;
      break;

    case SWITCH_DOWN: // run
      selected_mode = RC_mode::LOCOMOTION; // locomotion by default


      // stand mode
      if(left_select == SWITCH_UP && right_select == SWITCH_UP) {
        selected_mode = RC_mode::QP_STAND;
      }

      if(backflip_prep_edge_trigger.trigger(backflip_prepare_switch) && backflip_prepare_switch == SWITCH_MIDDLE) {
        initial_backflip_go_switch = backflip_go_switch;
      }

      // backflip
      if(backflip_prepare_switch == SWITCH_MIDDLE) {
        selected_mode = RC_mode::BACKFLIP_PRE;

        if(backflip_go_switch == SWITCH_DOWN && initial_backflip_go_switch != SWITCH_DOWN) {
          selected_mode = RC_mode::BACKFLIP;
        } else if(backflip_go_switch == SWITCH_UP) {
          initial_backflip_go_switch = SWITCH_UP;
        }
      }

      // gait selection
      int gait_id = left_select * 3 + right_select;
      constexpr int gait_table[9] = {0, //stand
                                     0, // trot
                                     1, // bounding
                                     2, // pronking
                                     3, // gallop
                                     5, // trot run
                                     6, // walk};
                                     7, // walk2?
                                     8, // pace
      };

      if(selected_mode == RC_mode::LOCOMOTION) {
        main_control_settings.variable[0] = gait_table[gait_id];
        main_control_settings.v_des[0] = v_scale * data.left_stick[1] * 0.5;
        main_control_settings.v_des[1] = v_scale * data.left_stick[0] * -1.;
        main_control_settings.v_des[2] = 0;
        main_control_settings.p_des[2] = 0.27 + 0.08 * data.knobs[1]; // todo?

        main_control_settings.omega_des[0] = 0;
        main_control_settings.omega_des[1] = 0;
        main_control_settings.omega_des[2] = -v_scale * data.right_stick[0];

        main_control_settings.rpy_des[1] = v_scale * data.right_stick[0];
        main_control_settings.alexa_mode = data.right_upper_switch == SWITCH_DOWN;

      } else if(selected_mode == RC_mode::QP_STAND) {
        main_control_settings.rpy_des[0] = data.left_stick[0] * 1.4;
        main_control_settings.rpy_des[1] = data.right_stick[1] * 0.46;
        main_control_settings.rpy_des[2] = -data.right_stick[0];

        main_control_settings.p_des[0] = 0;
        main_control_settings.p_des[1] = -0.667 * main_control_settings.rpy_des[0];
        main_control_settings.p_des[2] = data.left_stick[1] * .12;
        main_control_settings.alexa_mode = data.right_upper_switch == SWITCH_DOWN;
      }

      break;
  }

  bool trigger = mode_edge_trigger.trigger(selected_mode);
  if(trigger || selected_mode == RC_mode::OFF || selected_mode == RC_mode::RECOVERY_STAND) {
    if(trigger) {
      printf("MODE TRIGGER!\n");
    }
    main_control_settings.mode = selected_mode;
  }

}

/**
 * @brief      Function which handles the completion of an SBUS Packet and
 * overrides the LCM control settings as desired.
 */
void sbus_packet_complete_old() {
  pthread_mutex_lock(&lcm_get_set_mutex);

  int ch1 = read_sbus_channel(0);
  int ch2 = read_sbus_channel(1);
  int ch3 = read_sbus_channel(2);
  int ch4 = read_sbus_channel(3);

  // TODO: figure out what happens in channel 1
  int ch5 = read_sbus_channel(4);
  ch1 = ch5;

  int ch7 = read_sbus_channel(6);
  int ch8 = read_sbus_channel(7);
  int ch9 = read_sbus_channel(8);
  int ch10 = read_sbus_channel(9);
  int ch11 = read_sbus_channel(10);
  int ch12 = read_sbus_channel(11);
  int ch13 = read_sbus_channel(12);
  int ch15 = read_sbus_channel(14); // SE

  //printf("got sbus ch11: %d ch10 %d\n", ch11, ch10);

  //velocity scales between 1.0 and 3.0 (when 2.0f is used)
  // velocity scales between 1.0 and 1.5 (middle), 2.0 (right ends)
  float v_scale = 1.0f * (((float)(ch7 - 172)) / 1811.0f) + 1.0f;
  // Ignore commands if switched
  if (ch11 != 1811) {
    if (ch10 == 172) { // oh shit switch (OFF)
      main_control_settings.mode = RC_mode::OFF;

    } else if (ch10 == 992) { // Stand up recovery
      main_control_settings.mode = RC_mode::RECOVERY_STAND;

    } else if (ch10 == 1811) {  // ESTOP bar is down (Controller ON)
      main_control_settings.mode = RC_mode::LOCOMOTION; //Locomotion mode

      if(ch15==992){
        main_control_settings.mode = RC_mode::BACKFLIP_PRE; }
      else if(ch15==1811){
        main_control_settings.mode = RC_mode::BACKFLIP; }
      else if(ch12 == 1811){ 
        main_control_settings.mode = RC_mode::QP_STAND; }
      else if(ch12 == 992){
        main_control_settings.mode = RC_mode::VISION; }
    }

    // Use the joysticks for velocity and yaw control in locomotion gaits
    if (main_control_settings.mode == RC_mode::LOCOMOTION ||
        main_control_settings.mode == RC_mode::VISION) {
      // Analog channels return a value between ~200 and 1800, with sticks
      // centered at 1000

      if (ch9 == 1811) {
        main_control_settings.variable[0] = 4;
      }  // Stand
      else if (ch9 == 992 && (ch13==1811)) {
        main_control_settings.variable[0] = 0;
      }  // Trot
      else if (ch9 == 992 && (ch13==992)) {
        main_control_settings.variable[0] = 3;
      }  // Gallop
      else if (ch9 == 992 && (ch13==172)) {
        main_control_settings.variable[0] = 8;
      }  // Pacing 
      else if ((ch9 == 172) && (ch13 == 1811)) {
        main_control_settings.variable[0] = 5;
      }  // Trot run
      else if ((ch9 == 172) && (ch13 == 992)) {
        main_control_settings.variable[0] = 2;
      }  // Pronk
      else if ((ch9 == 172) && (ch13 == 172)) {
        main_control_settings.variable[0] = 1;
      }  // Bounding
      //main_control_settings.rpy_des[0] = ((float)ch4 - 1000) * .001f;
      //main_control_settings.rpy_des[1] = v_scale * ((float)ch1 - 1000) * .001f;
      //main_control_settings.rpy_des[2] = ((float)ch2 - 1000) * .001f;
      //main_control_settings.p_des[2] =   ((float)ch3 - 1000) * .001f;

      //printf("ch1, 2, 3, 4, 5: %d, %d, %d, %d , %d\n", ch1, ch2, ch3, ch4, ch5);
      main_control_settings.v_des[0] = v_scale * ((float)ch1-1000)*.001f;
      main_control_settings.v_des[1] = -v_scale *((float)ch4-1000)*.001f;
      main_control_settings.v_des[2] = 0;
      main_control_settings.p_des[2] = 0.25 + ((float)ch8 - 1000)*.0001f;;
      //printf("v scale : %f \n", v_scale);
      //printf("v des: %f, %f \n", main_control_settings.v_des[0], main_control_settings.v_des[1]);
      //printf("v scaled des: %f, %f \n", main_control_settings.v_des[0], main_control_settings.v_des[1]);

      main_control_settings.omega_des[0] = 0;
      main_control_settings.omega_des[1] = 0;
      main_control_settings.omega_des[2] = -v_scale * ((float)ch2 - 1000)*.002f;

      //      printf("%.3f %.3f %.3f %.3f\n",
      //             main_control_settings.rpy_des[0],
      //             main_control_settings.rpy_des[1],
      //             main_control_settings.rpy_des[2],
      //             main_control_settings.p_des[2]);


      // If using the stand gait, control orientations
//      if (main_control_settings.variable[0] == 4) {
//        main_control_settings.rpy_des[0] = ((float)ch4 - 1000) * .001f;
//        main_control_settings.rpy_des[1] = ((float)ch3 - 1000) * .001f;
//        main_control_settings.rpy_des[2] = ((float)ch2 - 1000) * .001f;
//
//        main_control_settings.p_des[0] = 0;
//        main_control_settings.p_des[1] = 0;
//        main_control_settings.p_des[2] = 0.21 + ((float)ch1 - 1000) * .0001f;
//      }
//      // For all other gaits, control the velocities
//      else {
//        main_control_settings.v_des[0] = v_scale * ((float)ch1 - 1000) * .002f;
//        main_control_settings.v_des[1] =
//            -v_scale * ((float)ch4 - 1000) * .0005f;
//        main_control_settings.v_des[2] = 0;
//        main_control_settings.p_des[2] = 0.21;
//
//        main_control_settings.omega_des[0] = 0;
//        main_control_settings.omega_des[1] = 0;
//        main_control_settings.omega_des[2] =
//            -v_scale * ((float)ch2 - 1000) * .003f;
//      }
    }
    // For standing modes (mpc or balance qp) control orientations
    else if (main_control_settings.mode == RC_mode::QP_STAND) {
      //main_control_settings.rpy_des[0] = ((float)ch4 - 1000) * .001f;
      //main_control_settings.rpy_des[1] = ((float)ch3 - 1000) * .001f;
      //main_control_settings.rpy_des[2] = -((float)ch2 - 1000)*.0015f;

      //main_control_settings.p_des[0] = 0;
      //main_control_settings.p_des[1] = 0;
      //main_control_settings.p_des[2] = 0.21 + ((float)ch1 - 1000) * .0001f;

      main_control_settings.rpy_des[0] = ((float)ch4-1000)*.0015f;
      main_control_settings.rpy_des[1] = ((float)ch3 - 1000)*.0005f;
      main_control_settings.rpy_des[2] = -((float)ch2 - 1000)*.0015f;

      main_control_settings.p_des[0] = 0;
      main_control_settings.p_des[1] = ((float)ch4 - 1000)*.00007f;;
      main_control_settings.p_des[2] = 0.25 + ((float)ch8 - 1000)*.0001f;
    }

    // For two contact standing mode
    else if (main_control_settings.mode == RC_mode::BACKFLIP) {
      main_control_settings.rpy_des[0] = ((float)ch4-1000)*.0015f;
      main_control_settings.rpy_des[1] = ((float)ch3 - 1000)*.0005f;
      main_control_settings.rpy_des[2] = -((float)ch2 - 1000)*.0015f;

      main_control_settings.p_des[0] = v_scale * ((float)ch1-1000)*.001f;
      main_control_settings.p_des[1] = -v_scale *((float)ch4-1000)*.001f;
      main_control_settings.p_des[2] = 0.0;
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
