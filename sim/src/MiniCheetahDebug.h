#ifndef MINICHEETAHDEBUG_H
#define MINICHEETAHDEBUG_H

#include <QDialog>
#include <cppTypes.h>

struct MiniCheetahDebugCommand {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Vec3<float> kp[4];
  Vec3<float> kd[4];
  Vec3<float> qd[4];
  bool enable[4][3];
};

struct MiniCheetahDebugData {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Vec3<float> p[4];
  Vec3<float> v[4];
};

namespace Ui {
class MiniCheetahDebug;
}

class MiniCheetahDebug : public QDialog
{
    Q_OBJECT

public:
    explicit MiniCheetahDebug(QWidget *parent = nullptr);
    void getDebugCommand(MiniCheetahDebugCommand& cmd);
    bool setDebugData(MiniCheetahDebugData& data);
    ~MiniCheetahDebug();

private:
    Ui::MiniCheetahDebug *ui;

    u32 update_counter = 0;
};

#endif // MINICHEETAHDEBUG_H
