#include <pthread.h>
#include <rt/rt_rc_interface.h>
#include "Utilities/EdgeTrigger.h"
#include <string.h> // memcpy
#include <stdio.h>
#include <rt/rt_sbus.h>
static pthread_mutex_t lcm_get_set_mutex =
PTHREAD_MUTEX_INITIALIZER; /**< mutex to protect gui settings coming over
                             LCM */

// Controller Settings
rc_control_settings rc_control;

/* ------------------------- HANDLERS ------------------------- */

// Controller Settings
void get_rc_control_settings(void *settings) {
  pthread_mutex_lock(&lcm_get_set_mutex);
  v_memcpy(settings, &rc_control, sizeof(rc_control_settings));
  pthread_mutex_unlock(&lcm_get_set_mutex);
}

//void get_rc_channels(void *settings) {
//pthread_mutex_lock(&lcm_get_set_mutex);
//v_memcpy(settings, &rc_channels, sizeof(rc_channels));
//pthread_mutex_unlock(&lcm_get_set_mutex);
//}

EdgeTrigger<int> mode_edge_trigger(0);
EdgeTrigger<TaranisSwitchState> backflip_prep_edge_trigger(SWITCH_UP);
EdgeTrigger<TaranisSwitchState> experiment_prep_edge_trigger(SWITCH_UP);
TaranisSwitchState initial_mode_go_switch = SWITCH_DOWN;

void sbus_packet_complete() {
  Taranis_X7_data data;
  update_taranis_x7(&data);

  float v_scale = data.knobs[0]*1.5f + 2.0f; // from 0.5 to 3.5
  float w_scale = 2.*v_scale; // from 1.0 to 7.0
  //printf("v scale: %f\n", v_scale);

  auto estop_switch = data.right_lower_right_switch;
  auto mode_selection_switch = data.left_lower_left_switch;
  auto mode_go_switch = data.left_upper_switch;

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

      if(backflip_prep_edge_trigger.trigger(mode_selection_switch) 
          && mode_selection_switch == SWITCH_MIDDLE) {
        initial_mode_go_switch = mode_go_switch;
      }

      // Experiment mode (two leg stance, vision, ...)
      if(experiment_prep_edge_trigger.trigger(mode_selection_switch) 
          && mode_selection_switch == SWITCH_DOWN) {
        initial_mode_go_switch = mode_go_switch;
      }


      // backflip
      if(mode_selection_switch == SWITCH_MIDDLE) {
        selected_mode = RC_mode::BACKFLIP_PRE;

        if(mode_go_switch == SWITCH_DOWN && initial_mode_go_switch != SWITCH_DOWN) {
          selected_mode = RC_mode::BACKFLIP;
        } else if(mode_go_switch == SWITCH_UP) {
          initial_mode_go_switch = SWITCH_UP;
        }
      } // Experiment Mode
      else if(mode_selection_switch == SWITCH_DOWN){
        int mode_id = left_select * 3 + right_select;

        if(mode_id == 0){ // Two leg stance
          selected_mode = RC_mode::TWO_LEG_STANCE_PRE;
          if(mode_go_switch == SWITCH_DOWN && initial_mode_go_switch != SWITCH_DOWN) {
            selected_mode = RC_mode::TWO_LEG_STANCE;
          } else if(mode_go_switch == SWITCH_UP) {
            initial_mode_go_switch = SWITCH_UP;
          }
        }
        else if(mode_id == 1){ // Vision 
          selected_mode = RC_mode::VISION;
        }
      }

      // gait selection
      int mode_id = left_select * 3 + right_select;

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


  // Deadband
  for(int i(0); i<2; ++i){
    data.left_stick[i] = deadband(data.left_stick[i], 0.1, -1., 1.);
    data.right_stick[i] = deadband(data.right_stick[i], 0.1, -1., 1.);
  }

  if(selected_mode == RC_mode::LOCOMOTION 
      || selected_mode == RC_mode::VISION) {
    rc_control.variable[0] = gait_table[mode_id];
    //rc_control.v_des[0] = v_scale * data.left_stick[1] * 0.5;
    //rc_control.v_des[1] = v_scale * data.left_stick[0] * -1.;
    rc_control.v_des[0] = v_scale * data.left_stick[1];
    rc_control.v_des[1] = -v_scale * data.left_stick[0];
    rc_control.v_des[2] = 0;

    rc_control.height_variation = data.knobs[1];
    //rc_control.p_des[2] = 0.27 + 0.08 * data.knobs[1]; // todo?

    rc_control.omega_des[0] = 0;
    rc_control.omega_des[1] = 0;
    rc_control.omega_des[2] = w_scale * data.right_stick[0];
    //rc_control.omega_des[2] = -v_scale * data.right_stick[0];

  } else if(selected_mode == RC_mode::QP_STAND || 
      selected_mode == RC_mode::TWO_LEG_STANCE) {
    //rc_control.rpy_des[0] = data.left_stick[0] * 1.4;
    //rc_control.rpy_des[1] = data.right_stick[1] * 0.46;
    rc_control.rpy_des[0] = data.left_stick[0];
    rc_control.rpy_des[1] = data.right_stick[1];
    rc_control.rpy_des[2] = data.right_stick[0];

    rc_control.height_variation = data.left_stick[1];

    rc_control.omega_des[0] = 0;
    rc_control.omega_des[1] = 0;
    rc_control.omega_des[2] = 0;
    //rc_control.p_des[1] = -0.667 * rc_control.rpy_des[0];
    //rc_control.p_des[2] = data.left_stick[1] * .12;
  } 
  break;
}

bool trigger = mode_edge_trigger.trigger(selected_mode);
if(trigger || selected_mode == RC_mode::OFF || selected_mode == RC_mode::RECOVERY_STAND) {
  if(trigger) {
    printf("MODE TRIGGER!\n");
  }
  rc_control.mode = selected_mode;
}

}

void *v_memcpy(void *dest, volatile void *src, size_t n) {
  void *src_2 = (void *)src;
  return memcpy(dest, src_2, n);
}

float deadband(float command, float deadbandRegion, float minVal, float maxVal){
  if (command < deadbandRegion && command > -deadbandRegion) {
    return 0.0;
  } else {
    return (command / (2)) * (maxVal - minVal);
  }
}


