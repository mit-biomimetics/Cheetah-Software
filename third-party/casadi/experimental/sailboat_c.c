/*CASADIMETA
:hull_resistance_N_IN 1
:hull_resistance_N_OUT 1
:hull_resistance_NAME_IN[0] U
:hull_resistance_NAME_OUT[0] X_0
*/

/*CASADIEXTERNAL hull_resistance_simple inline=1 */
extern "C" void hull_resistance_simple(const real_t* arg, real_t* res) {
  //real_t U = arg[0];
  res[0] = 0; // Neglected
}
/*CASADIEXTERNAL */
