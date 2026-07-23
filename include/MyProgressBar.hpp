#ifndef MYPROGRESSBAR_H
#define MYPROGRESSBAR_H

#include <cstdint>
#include <chrono>
#include <iostream>

#ifndef PROGRESS_BAR_COLUMN
#define PROGRESS_BAR_COLUMN 40
#endif

class MyProgressBar
{
private:
    bool is_done = false;
    uint32_t ticks = 0;
    uint32_t total_ticks = 0;

public:
    MyProgressBar(uint32_t total) : total_ticks{total} {}

    uint32_t operator++() { return ++ticks; }
    MyProgressBar &operator+=(const uint32_t tick)
    {
        ticks += tick;
        return *this;
    }

    void update(const uint32_t tick) { ticks = tick; }

    void reset()
    {
        ticks = 0;
        is_done = false;
    }

    void display()
    {
        uint8_t num_bars = static_cast<uint8_t>(PROGRESS_BAR_COLUMN * ticks / total_ticks);
        std::cout << "[";
        for (uint8_t i = 0; i < PROGRESS_BAR_COLUMN; i++)
        {
            if (i < num_bars)
                std::cout << "=";
            else if (i == num_bars)
                std::cout << ">";
            else
                std::cout << " ";
        }
        std::cout << "]";
    }

    void done()
    {
        if (is_done)
            return;
        is_done = true;
        display();
        std::cout << std::endl;
    }
};

#endif // MYPROGRESSBAR_H