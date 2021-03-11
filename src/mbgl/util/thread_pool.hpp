#pragma once

#include <mbgl/actor/mailbox.hpp>
#include <mbgl/actor/scheduler.hpp>

#include <algorithm>
#include <array>
#include <condition_variable>
#include <mutex>
#include <queue>
#include <thread>

namespace mbgl {

class ThreadedSchedulerBase : public Scheduler {
public:
    void schedule(std::function<void()>) override;

protected:
    ThreadedSchedulerBase() = default;
    ~ThreadedSchedulerBase() override;

    void terminate();
    std::thread makeSchedulerThread(size_t index);

    std::queue<std::function<void()>> queue;
    std::mutex mutex;
    std::condition_variable cv;
    bool terminated{false};
};

const size_t MAX_BACKGROUND_THREADS = 64;

/**
 * @brief ThreadScheduler implements Scheduler interface using a lightweight event loop
 *
 * @param N number of threads
 *
 * Note: If N == 1 all scheduled tasks are guaranteed to execute consequently;
 * otherwise, some of the scheduled tasks might be executed in parallel.
 */
class ThreadedScheduler : public ThreadedSchedulerBase {
public:
    ThreadedScheduler(std::size_t desiredThreads) {
        nThreads = std::min(desiredThreads, MAX_BACKGROUND_THREADS);

        for (std::size_t i = 0u; i < nThreads; ++i) {
            threads[i] = makeSchedulerThread(i);
        }
    }

    ~ThreadedScheduler() override {
        terminate();

        for (std::size_t i = 0; i < nThreads; i++) {
            auto& thread = threads[i];
            assert(std::this_thread::get_id() != thread.get_id());
            thread.join();
        }
    }

    mapbox::base::WeakPtr<Scheduler> makeWeakPtr() override { return weakFactory.makeWeakPtr(); }

private:
    std::size_t nThreads;
    std::array<std::thread, MAX_BACKGROUND_THREADS> threads;
    mapbox::base::WeakPtrFactory<Scheduler> weakFactory{this};
};

class SequencedScheduler : public ThreadedScheduler {
public:
    SequencedScheduler() : ThreadedScheduler(1) {}
};

class ThreadPool : public ThreadedScheduler {
public:
    ThreadPool(std::size_t desiredThreads) : ThreadedScheduler(desiredThreads) {}

};

} // namespace mbgl
