#include <cstdio>
#include <thread>
#include <unistd.h>

#include "rt/rt_sbus.h"
#include "rt/rt_rc_interface.h"

static struct {
 double     mode;
 double     impedance_scale;
 double     enable;
 double     emergency_damp;
 double     zero_leg[4];
 double     p_des[3];
 double     v_des[3];
 double     rpy_des[3];
 double     omega_des[3];
 double     p_des_slew_min[3];
 double     p_des_slew_max[3];
 double     rpy_des_slew_max[3];
 double     v_des_slew_min[3];
 double     v_des_slew_max[3];
 double     omegab_des_slew_max[3];
 double     emergency_damp_kd;
 double     alexa_mode;
 double     rc_configured;
 double     bonus_knee_torque;
 double     variable[3];
 double     want_cheater_mode;
} main_control_settings;

const char* switch_names[3] = {
    "UP", "--", "DN"
};

/*
     [UP]                        [UP]
  [UP]  [UP]  [ 0.1]  [-1.0]  [UP]  [UP]
         [-0.0 -0.0]  [ 0.0  0.0]
 */
void print_sbus(Taranis_X7_data& data) {
  printf("\n\n\n"
         "     [%s]                        [%s]\n"
         "  [%s]  [%s]  [%4.1f]  [%4.1f]  [%s]  [%s]\n"
         "         [%4.1f %4.1f]  [%4.1f %4.1f]\n",
         switch_names[(uint32_t)data.left_upper_switch],
         switch_names[(uint32_t)data.right_upper_switch],
         switch_names[(uint32_t)data.left_lower_left_switch],
         switch_names[(uint32_t)data.left_lower_right_switch],
         data.knobs[0], data.knobs[1],
         switch_names[(uint32_t)data.right_lower_left_switch],
         switch_names[(uint32_t)data.right_lower_right_switch],
         data.left_stick[0], data.left_stick[1], data.right_stick[0], data.right_stick[1]
  );
}

void print_left_y() {
 // printf("[%6.3f]\n", data.left_stick[1]);
  printf("cmd: %.3f\n", main_control_settings.v_des[0]);
  int cmd_int = (main_control_settings.v_des[0] * 50) + 50;
  for(int i = 0; i < cmd_int; i++) {
    printf("*");
  }
  printf("\n");
  printf("mode: %d\n\n", (int)main_control_settings.mode);
}

void sbus_packet_complete() {
  Taranis_X7_data data;
  update_taranis_x7(&data);

  float v_scale = data.knobs[0]*1.2f + 1.2f; // from 0.5 to 1.5

  auto estop_switch = data.right_lower_right_switch;
  auto backflip_prepare_switch = data.left_lower_left_switch;
  auto backflip_go_switch = data.left_upper_switch;

  auto left_select = data.left_lower_right_switch;
  auto right_select = data.right_lower_left_switch;

  switch(estop_switch) {

    case SWITCH_UP: // ESTOP
      main_control_settings.mode = RC_mode::OFF;
      break;

    case SWITCH_MIDDLE: // recover
      main_control_settings.mode = RC_mode::RECOVERY_STAND;
      break;

    case SWITCH_DOWN: // run
      main_control_settings.mode = RC_mode::LOCOMOTION; // locomotion by default


      // stand mode
      if(left_select == SWITCH_UP && right_select == SWITCH_UP) {
        main_control_settings.mode = RC_mode::QP_STAND;
      }

      // backflip
      if(backflip_prepare_switch == SWITCH_MIDDLE) {
        main_control_settings.mode = RC_mode::BACKFLIP_PRE;

        if(backflip_go_switch == SWITCH_DOWN) {
          main_control_settings.mode = RC_mode::BACKFLIP;
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

      if(main_control_settings.mode == RC_mode::LOCOMOTION) {
        main_control_settings.variable[0] = gait_table[gait_id];
        main_control_settings.v_des[0] = v_scale * data.left_stick[1] * 0.5;
        main_control_settings.v_des[1] = v_scale * data.left_stick[0] * -1.;
        main_control_settings.v_des[2] = 0;
        main_control_settings.p_des[2] = 0.25; // todo?

        main_control_settings.omega_des[0] = 0;
        main_control_settings.omega_des[1] = 0;
        main_control_settings.omega_des[2] = -v_scale * data.right_stick[0];

        main_control_settings.rpy_des[1] = v_scale * data.right_stick[0];

      } else if(main_control_settings.mode == RC_mode::QP_STAND) {
        main_control_settings.rpy_des[0] = data.left_stick[0] * 1.4;
        main_control_settings.rpy_des[1] = data.right_stick[1] * 0.46;
        main_control_settings.rpy_des[2] = -data.right_stick[0];

        main_control_settings.p_des[0] = 0;
        main_control_settings.p_des[1] = -0.667 * main_control_settings.rpy_des[0];
        main_control_settings.p_des[2] = data.left_stick[1] * .12;
      }




      break;
  }

}

#include <cmath>
static double last_cmd = 0.;
void check_d_rate() {
  if(std::abs(last_cmd - main_control_settings.v_des[0]) > 0.2) {
    printf("RATE error! (%6.3f -> %6.3f)\n", last_cmd, main_control_settings.v_des[0]);
  }
  last_cmd = main_control_settings.v_des[0];
}


std::thread lag_threads[12];

void lag() {
  for(auto& x : lag_threads) {
    x = std::thread([](){
      int64_t beans = 0;
      for(int64_t i = 2; i < 1000000000; i++) {
        beans += i * i - i + 1223123123/i + sqrt((double)i);
      }
      printf("done %ld\n", beans);
    });
  }
}

int main(int argc, char** argv) {
  (void)argc;
  (void)argv;

  int use_computer_port = 1;

  if(argc == 2 && argv[1][0] == 'r') {
    use_computer_port = 0;
  }

  lag();

  auto port = init_sbus(use_computer_port);
  // Taranis_X7_data controller_data;

  auto sbus_thread = std::thread([&](){
    while(true) {
      receive_sbus(port);
      if (receive_sbus(port)) {
        sbus_packet_complete();
        //update_taranis_x7(&controller_data);
      }
    }
  });

  for(;;) {
    for(int i = 0; i < 30; i++) {
      usleep(1000);
      check_d_rate();
    }
    print_left_y();
  }


  sbus_thread.join();


  return 0;
}
