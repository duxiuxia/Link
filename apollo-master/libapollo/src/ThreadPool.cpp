#include "ThreadPool.hpp"


template<>
ThreadPool<FIFO_POLICY>::ThreadPool (size_t numThreads){
    for(size_t i = 0 ; i < numThreads; i++){
        mWorkers.emplace_back(std::thread(
            [this] {
                while(true){
                    std::unique_lock<std::mutex> lock(this->queue_mutex);
                    while( this->isActive.load() && this->mTasks.empty())
                        this->condition.wait(lock);
                    if( ! this->isActive.load() && this->mTasks.empty())
                        return;
                    std::function<void(void)> lNextTask(this->mTasks.front());
                    this->mTasks.pop();
                    lock.unlock();
                    lNextTask();
                }
            }
        ));
    }
}
