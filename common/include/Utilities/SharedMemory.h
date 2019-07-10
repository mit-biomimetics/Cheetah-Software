/*! @file SharedMemory.h
 *  @brief Shared memory utilities for connecting the simulator program to the
 * robot program
 *
 *
 */
#ifndef PROJECT_SHAREDMEMORY_H
#define PROJECT_SHAREDMEMORY_H

#include <fcntl.h> /* For O_* constants */
#include <semaphore.h>
#include <sys/mman.h>
#include <sys/stat.h> /* For mode constants */
#include <unistd.h>
#include <cassert>
#include <cstring>
#include <stdexcept>
#include <string>
#include "cTypes.h"

#define DEVELOPMENT_SIMULATOR_SHARED_MEMORY_NAME "development-simulator"

/*!
 * A POSIX semaphore for shared memory.
 * See https://linux.die.net/man/7/sem_overview for more deatils
 */
class SharedMemorySemaphore {
 public:
  /*!
   * If semaphore is unitialized, initialize it and set its value.  This can be
   * called as many times as you want safely. It must be called at least once.
   * Only one process needs to call this, even if it is used in multiple
   * processes.
   *
   * Note that if init() is called after the semaphore has been initialized, it
   * will not change its value.
   * @param value The initial value of the semaphore.
   */
  void init(unsigned int value) {
    if (!_init) {
      if (sem_init(&_sem, 1, value)) {
        printf("[ERROR] Failed to initialize shared memory semaphore: %s\n",
               strerror(errno));
      } else {
        _init = true;
      }
    }
  }

  /*!
   * Increment the value of the semaphore.
   */
  void increment() { sem_post(&_sem); }

  /*!
   * If the semaphore's value is > 0, decrement the value.
   * Otherwise, wait until its value is > 0, then decrement.
   */
  void decrement() { sem_wait(&_sem); }

  /*!
   * If the semaphore's value is > 0, decrement the value and return true
   * Otherwise, return false (doesn't decrement or wait)
   * @return
   */
  bool tryDecrement() { return (sem_trywait(&_sem)) == 0; }

  /*!
   * Like decrement, but after waiting ms milliseconds, will give up
   * Returns true if the semaphore is successfully decremented
   */
  bool decrementTimeout(u64 seconds, u64 nanoseconds) {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    ts.tv_nsec += nanoseconds;
    ts.tv_sec += seconds;
    ts.tv_sec += ts.tv_nsec / 1000000000;
    ts.tv_nsec %= 1000000000;
    return (sem_timedwait(&_sem, &ts) == 0);
  }

  /*!
   * Delete the semaphore.  Note that deleting a semaphore in one process while
   * another is still using it results in very strange behavior.
   */
  void destroy() { sem_destroy(&_sem); }

 private:
  sem_t _sem;
  bool _init = false;
};

/*!
 * A container class for an object which is stored in shared memory.  This
 * object can then be viewed in multiple processes or programs.  Note that there
 * is significant overhead when creating a shared memory object, so it is
 * recommended that two programs that communicate should have one single large
 * SharedMemoryObject instead of many small ones.
 *
 * A name string is used to identify shared objects across different programs
 *
 * Before a shared memory object can be used, you must either allocate new
 * memory, or connect it to an existing shared memory object.
 *
 * Creating/deleting the memory can be done with createNew/closeNew.
 * Viewing an existing object allocated with createNew can be done with
 * attach/detach
 *
 * For an example, see test_sharedMemory.cpp
 */
template <typename T>
class SharedMemoryObject {
 public:
  SharedMemoryObject() = default;

  /*!
   * Allocate memory for the shared memory object and attach to it.
   * If allowOverwrite is true, and there's already an object with this name,
   * the old object is overwritten Note that if this happens, the object may be
   * initialized in a very weird state.
   *
   * Otherwise, if an object with the name already exists, throws a
   * std::runtime_error
   */
  bool createNew(const std::string& name, bool allowOverwrite = false) {
    bool hadToDelete = false;
    assert(!_data);
    _name = name;
    _size = sizeof(T);
    printf("[Shared Memory] open new %s, size %ld bytes\n", name.c_str(),
           _size);

    _fd = shm_open(name.c_str(), O_RDWR | O_CREAT,
                   S_IWUSR | S_IRUSR | S_IWGRP | S_IRGRP | S_IROTH);
    if (_fd == -1) {
      printf("[ERROR] SharedMemoryObject shm_open failed: %s\n",
             strerror(errno));
      throw std::runtime_error("Failed to create shared memory!");
      return false;
    }

    struct stat s;
    if (fstat(_fd, &s)) {
      printf("[ERROR] SharedMemoryObject::createNew(%s) stat: %s\n",
             name.c_str(), strerror(errno));
      throw std::runtime_error("Failed to create shared memory!");
      return false;
    }

    if (s.st_size) {
      printf(
          "[Shared Memory] SharedMemoryObject::createNew(%s) on something that "
          "wasn't new (size is %ld bytes)\n",
          _name.c_str(), s.st_size);
      hadToDelete = true;
      if (!allowOverwrite)
        throw std::runtime_error(
            "Failed to create shared memory - it already exists.");
      printf("\tusing existing shared memory!\n");
      // return false;
    }

    if (ftruncate(_fd, _size)) {
      printf("[ERROR] SharedMemoryObject::createNew(%s) ftruncate(%ld): %s\n",
             name.c_str(), _size, strerror(errno));
      throw std::runtime_error("Failed to create shared memory!");
      return false;
    }

    void* mem =
        mmap(nullptr, _size, PROT_READ | PROT_WRITE, MAP_SHARED, _fd, 0);
    if (mem == MAP_FAILED) {
      printf("[ERROR] SharedMemory::createNew(%s) mmap fail: %s\n",
             _name.c_str(), strerror(errno));
      throw std::runtime_error("Failed to create shared memory!");
      return false;
    }

    // there is a chance that the shared memory is not zeroed if we are reusing
    // old memory. this causes all sorts of weird issues, especially if the
    // layout of the object in memory has changed.
    memset(mem, 0, _size);

    _data = (T*)mem;
    return hadToDelete;
  }

  /*!
   * Attach to an existing shared memory object.
   */
  void attach(const std::string& name) {
    assert(!_data);
    _name = name;
    _size = sizeof(T);
    printf("[Shared Memory] open existing %s size %ld bytes\n", name.c_str(),
           _size);
    _fd = shm_open(name.c_str(), O_RDWR,
                   S_IWUSR | S_IRUSR | S_IWGRP | S_IRGRP | S_IROTH);
    if (_fd == -1) {
      printf("[ERROR] SharedMemoryObject::attach shm_open(%s) failed: %s\n",
             _name.c_str(), strerror(errno));
      throw std::runtime_error("Failed to create shared memory!");
      return;
    }

    struct stat s;
    if (fstat(_fd, &s)) {
      printf("[ERROR] SharedMemoryObject::attach(%s) stat: %s\n", name.c_str(),
             strerror(errno));
      throw std::runtime_error("Failed to create shared memory!");
      return;
    }

    if ((size_t)s.st_size != _size) {
      printf(
          "[ERROR] SharedMemoryObject::attach(%s) on something that was "
          "incorrectly "
          "sized (size is %ld bytes, should be %ld)\n",
          _name.c_str(), s.st_size, _size);
      throw std::runtime_error("Failed to create shared memory!");
      return;
    }

    void* mem =
        mmap(nullptr, _size, PROT_READ | PROT_WRITE, MAP_SHARED, _fd, 0);
    if (mem == MAP_FAILED) {
      printf("[ERROR] SharedMemory::attach(%s) mmap fail: %s\n", _name.c_str(),
             strerror(errno));
      throw std::runtime_error("Failed to create shared memory!");
      return;
    }

    _data = (T*)mem;
  }

  /*!
   * Free memory associated with the current open shared memory object.  The
   * object could have been opened with either attach or createNew.  After
   * calling this, no process can use this shared object
   */
  void closeNew() {
    assert(_data);
    // first, unmap
    if (munmap((void*)_data, _size)) {
      printf("[ERROR] SharedMemoryObject::closeNew (%s) munmap %s\n",
             _name.c_str(), strerror(errno));
      throw std::runtime_error("Failed to create shared memory!");
      return;
    }

    _data = nullptr;

    if (shm_unlink(_name.c_str())) {
      printf("[ERROR] SharedMemoryObject::closeNew (%s) shm_unlink %s\n",
             _name.c_str(), strerror(errno));
      throw std::runtime_error("Failed to create shared memory!");
      return;
    }

    // close fd
    if (close(_fd)) {
      printf("[ERROR] SharedMemoryObject::closeNew (%s) close %s\n",
             _name.c_str(), strerror(errno));
      throw std::runtime_error("Failed to create shared memory!");
      return;
    }

    _fd = 0;
  }

  /*!
   * Close this view of the currently opened shared memory object. The object
   * can be opened with either attach or createNew.  After calling this, this
   * process can no longer use this shared object, but other processes still
   * can.
   */
  void detach() {
    assert(_data);
    // first, unmap
    if (munmap((void*)_data, _size)) {
      printf("[ERROR] SharedMemoryObject::detach (%s) munmap %s\n",
             _name.c_str(), strerror(errno));
      throw std::runtime_error("Failed to create shared memory!");
      return;
    }

    _data = nullptr;

    // close fd
    if (close(_fd)) {
      printf("[ERROR] SharedMemoryObject::detach (%s) close %s\n",
             _name.c_str(), strerror(errno));
      throw std::runtime_error("Failed to create shared memory!");
      return;
    }

    _fd = 0;
  }

  /*!
   * Get the shared memory object.
   */
  T* get() {
    assert(_data);
    return _data;
  }

  /*!
   * Get the shared memory object.
   */
  T& operator()() {
    assert(_data);
    return *_data;
  }

 private:
  T* _data = nullptr;
  std::string _name;
  size_t _size;
  int _fd;
};

#endif  // PROJECT_SHAREDMEMORY_H
