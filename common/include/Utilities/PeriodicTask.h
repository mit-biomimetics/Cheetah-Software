/*!
 * @file PeriodicTask.h
 * @brief Implementation of a periodic function running in a separate thread.
 * Periodic tasks have a task manager, which measure how long they take to run.
 */

#ifndef PROJECT_PERIODICTASK_H
#define PROJECT_PERIODICTASK_H

#include <string>
#include <thread>
#include <vector>

class PeriodicTaskManager;

/*!
 * A single periodic task which will call run() at the given frequency
 */
class PeriodicTask {
 public:
  PeriodicTask(PeriodicTaskManager* taskManager, float period,
               std::string name);
  void start();
  void stop();
  void printStatus();
  void clearMax();
  bool isSlow();
  virtual void init() = 0;
  virtual void run() = 0;
  virtual void cleanup() = 0;
  virtual ~PeriodicTask() { stop(); }

  /*!
   * Get the desired period for the task
   */
  float getPeriod() { return _period; }

  /*!
   * Get how long the most recent run took
   */
  float getRuntime() { return _lastRuntime; }

  /*!
   * Get the maximum time in between runs
   */
  float getMaxPeriod() { return _maxPeriod; }

  /*!
   * Get the maximum time it took for a run
   */
  float getMaxRuntime() { return _maxRuntime; }

 private:
  void loopFunction();

  float _period;
  volatile bool _running = false;
  float _lastRuntime = 0;
  float _lastPeriodTime = 0;
  float _maxPeriod = 0;
  float _maxRuntime = 0;
  std::string _name;
  std::thread _thread;
};

/*!
 * A collection of periodic tasks which can be monitored together
 */
class PeriodicTaskManager {
 public:
  PeriodicTaskManager() = default;
  ~PeriodicTaskManager();
  void addTask(PeriodicTask* task);
  void printStatus();
  void printStatusOfSlowTasks();
  void stopAll();

 private:
  std::vector<PeriodicTask*> _tasks;
};

/*!
 * A periodic task for calling a function
 */
class PeriodicFunction : public PeriodicTask {
 public:
  PeriodicFunction(PeriodicTaskManager* taskManager, float period,
                   std::string name, void (*function)())
      : PeriodicTask(taskManager, period, name), _function(function) {}
  void cleanup() {}
  void init() {}
  void run() { _function(); }

  ~PeriodicFunction() { stop(); }

 private:
  void (*_function)() = nullptr;
};

/*!
 * A periodic task for printing the status of all tasks in the task manager
 */
class PrintTaskStatus : public PeriodicTask {
 public:
  PrintTaskStatus(PeriodicTaskManager* tm, float period)
      : PeriodicTask(tm, period, "print-tasks"), _tm(tm) {}
  void run() override { 
    // DH: Disable printing
    //_tm->printStatus(); 
  }

  void init() override {}

  void cleanup() override {}

 private:
  PeriodicTaskManager* _tm;
};

/*!
 * A periodic task for calling a member function
 */
template <typename T>
class PeriodicMemberFunction : public PeriodicTask {
 public:
  PeriodicMemberFunction(PeriodicTaskManager* taskManager, float period,
                         std::string name, void (T::*function)(), T* obj)
      : PeriodicTask(taskManager, period, name),
        _function(function),
        _obj(obj) {}

  void cleanup() {}
  void init() {}
  void run() { (_obj->*_function)(); }

 private:
  void (T::*_function)();
  T* _obj;
};

#endif  // PROJECT_PERIODICTASK_H
