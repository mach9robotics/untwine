/*****************************************************************************
 *   Copyright (c) 2026, Mach9 Robotics                                      *
 *                                                                           *
 *   Unit tests for untwine::ThreadPool -- the drain and error-propagation   *
 *   guarantees that bu::ChunkWriter's deferred compression relies on:       *
 *   every enqueued chunk must be written before stop() returns, the bounded *
 *   queue must apply backpressure without dropping work, and a worker       *
 *   failure must surface so stop() can rethrow it before the file is        *
 *   finalized.                                                              *
 ****************************************************************************/

#include <gtest/gtest.h>

#include <atomic>
#include <stdexcept>
#include <string>

#include "untwine/ThreadPool.hpp"

using namespace untwine;

// join() must run every enqueued task before returning. ChunkWriter::stop() depends on this so
// that no chunk is left uncompressed/unwritten when the COPC file is finalized.
TEST(ThreadPool, join_drains_all_tasks)
{
    ThreadPool pool(4);
    std::atomic<int> n{0};
    for (int i = 0; i < 10000; ++i)
        pool.add([&n]() { n.fetch_add(1, std::memory_order_relaxed); });
    pool.join();
    EXPECT_EQ(n.load(), 10000);
}

// A bounded queue (ChunkWriter uses depth 16) must apply backpressure -- add() blocks when full --
// without dropping tasks: every task still runs.
TEST(ThreadPool, bounded_queue_runs_every_task)
{
    ThreadPool pool(2, /*queueSize=*/4);
    std::atomic<int> n{0};
    for (int i = 0; i < 2000; ++i)
        pool.add([&n]() { n.fetch_add(1, std::memory_order_relaxed); });
    pool.join();
    EXPECT_EQ(n.load(), 2000);
}

// With trapping on, a throwing task is caught and surfaced via clearErrors, and the other tasks
// still run. This is how ChunkWriter::stop() turns a worker compression/write failure into a
// thrown FatalError.
TEST(ThreadPool, traps_and_surfaces_worker_errors)
{
    ThreadPool pool(2);
    pool.trap(true, "catchall");
    std::atomic<int> n{0};
    pool.add([]() { throw std::runtime_error("boom"); });
    pool.add([&n]() { n.fetch_add(1, std::memory_order_relaxed); });
    pool.join();

    std::vector<std::string> errs = pool.clearErrors();
    ASSERT_FALSE(errs.empty());
    bool found = false;
    for (const std::string& e : errs)
        if (e.find("boom") != std::string::npos)
            found = true;
    EXPECT_TRUE(found);
    EXPECT_EQ(n.load(), 1);  // the non-throwing task still ran
}

// A non-std::exception throw is caught and reported with the configured catchall message.
TEST(ThreadPool, traps_nonstd_exception_with_catchall)
{
    ThreadPool pool(2);
    pool.trap(true, "catchall-msg");
    pool.add([]() { throw 42; });
    pool.join();

    std::vector<std::string> errs = pool.clearErrors();
    ASSERT_FALSE(errs.empty());
    EXPECT_NE(errs[0].find("catchall-msg"), std::string::npos);
}
