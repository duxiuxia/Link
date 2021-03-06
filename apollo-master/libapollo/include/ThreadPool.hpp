/*
 *Copyright (c) 2012 Jakob Progsch
 *
 *This software is provided 'as-is', without any express or implied
 *warranty. In no event will the authors be held liable for any damages
 *arising from the use of this software.
 *
 *Permission is granted to anyone to use this software for any purpose,
 *including commercial applications, and to alter it and redistribute it
 *freely, subject to the following restrictions:
 *
 *  1. The origin of this software must not be misrepresented; you must not
 *  claim that you wrote the original software. If you use this software
 *  in a product, an acknowledgment in the product documentation would be
 *  appreciated but is not required.
 *
 *  2. Altered source versions must be plainly marked as such, and must not be
 *  misrepresented as being the original software.

 *  3. This notice may not be removed or altered from any source
 *  distribution.
 *
 * Copied from https://github.com/en4bz/ThreadPool
 */

#ifndef THREAD_POOL_HPP
#define THREAD_POOL_HPP

#include <stack>
#include <queue>
#include <mutex>
#include <memory>
#include <thread>
#include <vector>
#include <future>
#include <stdexcept>
#include <functional>
#include <type_traits>
#include <condition_variable>

#include <iostream>

class prioritized_task{
private:
    int priority;
    std::function<void(void)> task;
public:
    prioritized_task(int p, std::function<void(void)>&& f) : priority(p), task(f) {}

    bool operator< (const prioritized_task& other) const{
        return this->priority < other.priority;
    }

    void operator()(void){
        task();
    }
};

typedef std::queue<std::function<void(void)>>   FIFO_POLICY;
typedef std::stack<std::function<void(void)>>   LIFO_POLICY;
typedef std::priority_queue<prioritized_task>   PRIORITY_POLICY;

template <typename policy_type = FIFO_POLICY>
class ThreadPool {
private:
    std::vector<std::thread> mWorkers;
    policy_type mTasks;

    std::mutex queue_mutex;
    std::condition_variable condition;
    std::atomic<bool> isActive;
public:
    ThreadPool (size_t numThreads = std::thread::hardware_concurrency() ){
        for(size_t i = 0 ; i < numThreads; i++){
            mWorkers.emplace_back(std::thread(
                [this] {
                    while(true){
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        while( this->isActive.load() && this->mTasks.empty())
                            this->condition.wait(lock);
                        if( ! this->isActive.load() && this->mTasks.empty())
                            return;
                        auto lNextTask = std::move (this->mTasks.top());
                        this->mTasks.pop();
                        lock.unlock();
                        lNextTask();
                    }
                }
            ));
        }
    }

    ThreadPool(const ThreadPool& to_copy) = delete; //Probably wouldn't and shouldn't copy this.

    template<class F, class... Args> //Below is the return type...Yes it is ridiculous, but it works. Enabled if policy_tpye IS NOT PRIORITY
    typename std::enable_if< ! std::is_same<policy_type,PRIORITY_POLICY>::value, std::future<typename std::result_of<F(Args...)>::type> >::type
    enqueue(F&& f, Args&&... args){
        typedef typename std::result_of<F(Args...)>::type return_type;
        // Don't allow enqueueing after stopping the pool
        if ( ! isActive.load() )
            throw std::runtime_error("enqueue on stopped ThreadPool");

        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...) );

        std::future<return_type> result = task->get_future();
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            mTasks.push([task](){ (*task)(); });
        }
        condition.notify_one();
        return result;
    }

    template<class F, class... Args> //Below is the return type...Yes it is ridiculous, but it works. Enabled if policy_tpye IS PRIORITY
    typename std::enable_if< std::is_same<policy_type,PRIORITY_POLICY>::value, std::future<typename std::result_of<F(Args...)>::type> >::type
    enqueue(int priority, F&& f, Args&&... args){
        typedef typename std::result_of<F(Args...)>::type return_type;
        // Don't allow enqueueing after stopping the pool
        if ( ! isActive.load() )
            throw std::runtime_error("enqueue on stopped ThreadPool");

        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...) );

        std::future<return_type> res = task->get_future();
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            mTasks.push(prioritized_task(priority, [task](void){ (*task)();}));
        }
        condition.notify_one();
        return res;
    }

    int pending(void){
        std::unique_lock<std::mutex> lock(queue_mutex);
        return this->mTasks.size();
    }

    ~ThreadPool(void){
        this->isActive = false;
        condition.notify_all();
        for(std::thread& t : mWorkers)
        t.join();
    }
};

template<>
ThreadPool<FIFO_POLICY>::ThreadPool (size_t numThreads);
#endif
