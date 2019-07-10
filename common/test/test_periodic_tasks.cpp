#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Utilities/PeriodicTask.h"

class TestPeriodicTask : public PeriodicTask {
 public:
  using PeriodicTask::PeriodicTask;
  int _counter = 0;
  bool _cleanedUp = false;
  bool _init = false;
  bool _slow = false;

  void run() override {
    _counter++;
    if (_slow) usleep(15000);
  }

  void init() override { _init = true; }

  void cleanup() override { _cleanedUp = true; }
};

static int scount = 0;

void inc_scount() { scount++; }

TEST(PeriodicTask, test1) {
  PeriodicTaskManager taskManager;
  TestPeriodicTask task1(&taskManager, 0.01f, "test-task-1");
  TestPeriodicTask task2(&taskManager, 0.02f, "test-task-2");

  task1.start();
  task2.start();

  for (int i = 0; i < 10; i++) {
    // taskManager.printStatus();
    usleep(50000);
  }

  task1.start();
  task1.start();

  EXPECT_TRUE(task2._counter > 0);
  EXPECT_TRUE(task2._init);
  EXPECT_FALSE(task2._cleanedUp);

  EXPECT_TRUE(task2._counter > 0);
  EXPECT_TRUE(task2._init);
  EXPECT_FALSE(task2._cleanedUp);

  taskManager.stopAll();

  TestPeriodicTask task3(&taskManager, .01f, "slow-task");
  task3._slow = true;
  task3.start();

  for (int i = 0; i < 10; i++) {
    // taskManager.printStatus();
    usleep(50000);
  }

  taskManager.stopAll();
  taskManager.stopAll();

  PeriodicFunction pf(&taskManager, 0.01f, "func-test", &inc_scount);
  pf.start();

  usleep(100000);
  EXPECT_TRUE(scount > 0);
  printf("scount = %d\n", scount);

  // taskManager.stopAll(); test destructors cleaning things up instead
}
