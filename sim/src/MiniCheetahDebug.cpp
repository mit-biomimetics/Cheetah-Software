#include "MiniCheetahDebug.h"
#include "ui_MiniCheetahDebug.h"
//#include "../../build/sim/sim_autogen/include/ui_MiniCheetahDebug.h" // todo, remove me

constexpr float KP_MAX = 10;
constexpr float KD_MAX = 0.5;

constexpr float SLIDER_MIN = -1;
constexpr float SLIDER_MAX = 1;

static float coerce(float in, float min, float max) {
  if(in > max) return max;
  if(in < min) return min;
  return in;
}

static float slider_map(int in) {
  float span = SLIDER_MAX - SLIDER_MIN;
  return span * (float)(in) / 100.f - (span / 2);
}

static QString make_string(float v) {
  char buff[128];
  sprintf(buff, "%0.3f", v);
  return QString(buff);
}

MiniCheetahDebug::MiniCheetahDebug(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::MiniCheetahDebug)
{
    ui->setupUi(this);
    setWindowTitle("Mini Cheetah Debug Window");

    // hack
    MiniCheetahDebugData zero_data;
    for(auto& p : zero_data.p) p.setZero();
    for(auto& v : zero_data.v) v.setZero();
    setDebugData(zero_data);

  ui->abad_1_kp->setText("0");
  ui->hip_1_kp->setText("0");
  ui->knee_1_kp->setText("0");
  ui->abad_2_kp->setText("0");
  ui->hip_2_kp->setText("0");
  ui->knee_2_kp->setText("0");
  ui->abad_3_kp->setText("0");
  ui->hip_3_kp->setText("0");
  ui->knee_3_kp->setText("0");
  ui->abad_4_kp->setText("0");
  ui->hip_4_kp->setText("0");
  ui->knee_4_kp->setText("0");
  ui->abad_1_kd->setText("0");
  ui->hip_1_kd->setText("0");
  ui->knee_1_kd->setText("0");
  ui->abad_2_kd->setText("0");
  ui->hip_2_kd->setText("0");
  ui->knee_2_kd->setText("0");
  ui->abad_3_kd->setText("0");
  ui->hip_3_kd->setText("0");
  ui->knee_3_kd->setText("0");
  ui->abad_4_kd->setText("0");
  ui->hip_4_kd->setText("0");
  ui->knee_4_kd->setText("0");

  ui->abad_1_qd->setValue(50);
  ui->hip_1_qd->setValue(50);
  ui->knee_1_qd->setValue(50);
  ui->abad_2_qd->setValue(50);
  ui->hip_2_qd->setValue(50);
  ui->knee_2_qd->setValue(50);
  ui->abad_3_qd->setValue(50);
  ui->hip_3_qd->setValue(50);
  ui->knee_3_qd->setValue(50);
  ui->abad_4_qd->setValue(50);
  ui->hip_4_qd->setValue(50);
  ui->knee_4_qd->setValue(50);
}

MiniCheetahDebug::~MiniCheetahDebug()
{
    delete ui;
}

void MiniCheetahDebug::getDebugCommand(MiniCheetahDebugCommand &cmd) {
  cmd.enable[0][0] = ui->abad_1_enable->isChecked();
  cmd.enable[0][1] = ui->hip_1_enable->isChecked();
  cmd.enable[0][2] = ui->knee_1_enable->isChecked();

  cmd.enable[1][0] = ui->abad_2_enable->isChecked();
  cmd.enable[1][1] = ui->hip_2_enable->isChecked();
  cmd.enable[1][2] = ui->knee_2_enable->isChecked();

  cmd.enable[2][0] = ui->abad_3_enable->isChecked();
  cmd.enable[2][1] = ui->hip_3_enable->isChecked();
  cmd.enable[2][2] = ui->knee_3_enable->isChecked();

  cmd.enable[3][0] = ui->abad_4_enable->isChecked();
  cmd.enable[3][1] = ui->hip_4_enable->isChecked();
  cmd.enable[3][2] = ui->knee_4_enable->isChecked();

  cmd.kp[0][0] = coerce(ui->abad_1_kp->text().toFloat(), -KP_MAX, KP_MAX);
  cmd.kp[0][1] = coerce(ui->hip_1_kp->text().toFloat(), -KP_MAX, KP_MAX);
  cmd.kp[0][2] = coerce(ui->knee_1_kp->text().toFloat(), -KP_MAX, KP_MAX);

  cmd.kp[1][0] = coerce(ui->abad_2_kp->text().toFloat(), -KP_MAX, KP_MAX);
  cmd.kp[1][1] = coerce(ui->hip_2_kp->text().toFloat(), -KP_MAX, KP_MAX);
  cmd.kp[1][2] = coerce(ui->knee_2_kp->text().toFloat(), -KP_MAX, KP_MAX);

  cmd.kp[2][0] = coerce(ui->abad_3_kp->text().toFloat(), -KP_MAX, KP_MAX);
  cmd.kp[2][1] = coerce(ui->hip_3_kp->text().toFloat(), -KP_MAX, KP_MAX);
  cmd.kp[2][2] = coerce(ui->knee_3_kp->text().toFloat(), -KP_MAX, KP_MAX);

  cmd.kp[3][0] = coerce(ui->abad_4_kp->text().toFloat(), -KP_MAX, KP_MAX);
  cmd.kp[3][1] = coerce(ui->hip_4_kp->text().toFloat(), -KP_MAX, KP_MAX);
  cmd.kp[3][2] = coerce(ui->knee_4_kp->text().toFloat(), -KP_MAX, KP_MAX);

  cmd.kd[0][0] = coerce(ui->abad_1_kd->text().toFloat(), -KD_MAX, KD_MAX);
  cmd.kd[0][1] = coerce(ui->hip_1_kd->text().toFloat(), -KD_MAX, KD_MAX);
  cmd.kd[0][2] = coerce(ui->knee_1_kd->text().toFloat(), -KD_MAX, KD_MAX);

  cmd.kd[1][0] = coerce(ui->abad_2_kd->text().toFloat(), -KD_MAX, KD_MAX);
  cmd.kd[1][1] = coerce(ui->hip_2_kd->text().toFloat(), -KD_MAX, KD_MAX);
  cmd.kd[1][2] = coerce(ui->knee_2_kd->text().toFloat(), -KD_MAX, KD_MAX);

  cmd.kd[2][0] = coerce(ui->abad_3_kd->text().toFloat(), -KD_MAX, KD_MAX);
  cmd.kd[2][1] = coerce(ui->hip_3_kd->text().toFloat(), -KD_MAX, KD_MAX);
  cmd.kd[2][2] = coerce(ui->knee_3_kd->text().toFloat(), -KD_MAX, KD_MAX);

  cmd.kd[3][0] = coerce(ui->abad_4_kd->text().toFloat(), -KD_MAX, KD_MAX);
  cmd.kd[3][1] = coerce(ui->hip_4_kd->text().toFloat(), -KD_MAX, KD_MAX);
  cmd.kd[3][2] = coerce(ui->knee_4_kd->text().toFloat(), -KD_MAX, KD_MAX);

  cmd.qd[0][0] = slider_map(ui->abad_1_qd->value());
  cmd.qd[0][1] = slider_map(ui->hip_1_qd->value());
  cmd.qd[0][2] = slider_map(ui->knee_1_qd->value());

  cmd.qd[1][0] = slider_map(ui->abad_2_qd->value());
  cmd.qd[1][1] = slider_map(ui->hip_2_qd->value());
  cmd.qd[1][2] = slider_map(ui->knee_2_qd->value());

  cmd.qd[2][0] = slider_map(ui->abad_3_qd->value());
  cmd.qd[2][1] = slider_map(ui->hip_3_qd->value());
  cmd.qd[2][2] = slider_map(ui->knee_3_qd->value());

  cmd.qd[3][0] = slider_map(ui->abad_4_qd->value());
  cmd.qd[3][1] = slider_map(ui->hip_4_qd->value());
  cmd.qd[3][2] = slider_map(ui->knee_4_qd->value());
}

bool MiniCheetahDebug::setDebugData(MiniCheetahDebugData &msg) {
  update_counter++;
  if(update_counter >= 10) {
    update_counter = 0;
  } else {
    return false;
  }

  ui->abad_1_pos->setText(make_string(msg.p[0][0]));
  ui->hip_1_pos->setText(make_string(msg.p[0][1]));
  ui->knee_1_pos->setText(make_string(msg.p[0][2]));

  ui->abad_2_pos->setText(make_string(msg.p[1][0]));
  ui->hip_2_pos->setText(make_string(msg.p[1][1]));
  ui->knee_2_pos->setText(make_string(msg.p[1][2]));

  ui->abad_3_pos->setText(make_string(msg.p[2][0]));
  ui->hip_3_pos->setText(make_string(msg.p[2][1]));
  ui->knee_3_pos->setText(make_string(msg.p[2][2]));

  ui->abad_4_pos->setText(make_string(msg.p[3][0]));
  ui->hip_4_pos->setText(make_string(msg.p[3][1]));
  ui->knee_4_pos->setText(make_string(msg.p[3][2]));

  ui->abad_1_vel->setText(make_string(msg.v[0][0]));
  ui->hip_1_vel->setText(make_string(msg.v[0][1]));
  ui->knee_1_vel->setText(make_string(msg.v[0][2]));

  ui->abad_2_vel->setText(make_string(msg.v[1][0]));
  ui->hip_2_vel->setText(make_string(msg.v[1][1]));
  ui->knee_2_vel->setText(make_string(msg.v[1][2]));

  ui->abad_3_vel->setText(make_string(msg.v[2][0]));
  ui->hip_3_vel->setText(make_string(msg.v[2][1]));
  ui->knee_3_vel->setText(make_string(msg.v[2][2]));

  ui->abad_4_vel->setText(make_string(msg.v[3][0]));
  ui->hip_4_vel->setText(make_string(msg.v[3][1]));
  ui->knee_4_vel->setText(make_string(msg.v[3][2]));

  return true;
}