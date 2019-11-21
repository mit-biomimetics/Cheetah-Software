/*! @file test_SharedMemory.cpp
 *  @brief Test the shared memory
 *
 *
 */

#include "Utilities/SharedMemory.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace SharedMemoryTestNamespace {
struct TestType {
  u32 value;
  SharedMemorySemaphore semaphore;
};
}  // namespace SharedMemoryTestNamespace

using namespace SharedMemoryTestNamespace;

// test shared memory
TEST(SharedMemory, MemoryAndSemaphore) {
  SharedMemoryObject<TestType> sharedObject, oldStaleSharedMemory;
  std::string memoryName = "/shared-memory-test";

  // first open and close shared memory to reset everything
  // allow overwriting because the last run of the test could have left memory
  // open
  sharedObject.createNew(memoryName, true);
  EXPECT_TRUE(sharedObject.get());
  sharedObject.closeNew();

  // now we should be able to open shared memory without overwriting
  EXPECT_FALSE(oldStaleSharedMemory.createNew(memoryName, false));
  oldStaleSharedMemory.get()->value =
      12;  // this shouldn't segfault if the memory is mapped succesfully

  // now we should have to overwrite already existing shared memory
  EXPECT_TRUE(sharedObject.createNew(memoryName, true));
  sharedObject.closeNew();  // and close it

  // now we can open without overwriting
  EXPECT_FALSE(sharedObject.createNew(memoryName, false));

  const u32 valueFromParent = 2;
  const u32 valueFromChildSuccess = 3;
  const u32 valueFromChildFailure = 4;

  // value should be zeroed
  EXPECT_TRUE(0 == sharedObject.get()->value);
  // set up tosend value from parent to child process
  sharedObject.get()->value = valueFromParent;
  sharedObject.get()->semaphore.init(0);

  // create child process
  pid_t pid = fork();
  if (!pid) {
    // child runs this: open shared memory
    SharedMemoryObject<TestType> childView;
    childView.attach(memoryName);

    // check to see we got the correct value:
    if (childView().value == valueFromParent) {
      childView().value = valueFromChildSuccess;
    } else {
      // otherwise signal an error
      childView().value = valueFromChildFailure;
    }

    // inform the parent that it can now read the child value from shared memory
    childView().semaphore.increment();

    // close our view and exit child process
    childView.detach();
    exit(0);
  }

  sharedObject.get()->semaphore.decrement();  // wait until child writes
  EXPECT_TRUE(valueFromChildSuccess == sharedObject.get()->value);
  sharedObject.closeNew();
}