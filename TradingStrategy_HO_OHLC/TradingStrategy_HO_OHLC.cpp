// Compile with:
//   g++ -std=c++17 -O2 -o backtest main.cpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <ctime>
#include <limits>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <deque>

struct Bar {
    std::time_t t;
    double open, high, low, close, volume;
};

static std::time_t parse_datetime(const std::string& date_str, const std::string& time_str) {
    std::tm tm = {};
    int dd, mm, yyyy, HH, MM;
    char sep;

    {
        std::istringstream ss(date_str);
        ss >> mm >> sep >> dd >> sep >> yyyy;
        if (!ss || sep != '/') {
            std::cerr << "Date parse error: " << date_str << "\n";
            std::exit(1);
        }
    }

    {
        std::istringstream ss(time_str);
        ss >> HH >> sep >> MM;
        if (!ss || sep != ':') {
            std::cerr << "Time parse error: " << time_str << "\n";
            std::exit(1);
        }
    }

    tm.tm_year = yyyy - 1900;
    tm.tm_mon = mm - 1;
    tm.tm_mday = dd;
    tm.tm_hour = HH;
    tm.tm_min = MM;
    tm.tm_sec = 0;
    return std::mktime(&tm);
}

int main() {
    // Parameters
    const std::string CSV_PATH = "HO-5minHLV.csv";
    const int    BARS_BACK = 17001;
    const int    L = 12700;
    const double S = 0.010;
    const double SL_PG = 70.0;
    const double PV = 42000.0;
    const double E0 = 100000.0;
    const std::string IN_START_DATE = "01/01/1980";
    const std::string IN_END_DATE = "01/01/2000";
    const std::string OUT_END_DATE = "03/23/2023";

    // 1) Count lines
    std::ifstream file(CSV_PATH);
    if (!file.is_open()) {
        std::perror("Failed to open CSV");
        return 1;
    }
    std::string line;
    std::getline(file, line); // skip header
    int nBars = 0;
    while (std::getline(file, line)) ++nBars;
    file.clear();
    file.seekg(0);
    std::getline(file, line); // skip header again

    // 2) Allocate containers
    std::vector<Bar>    bars;    bars.reserve(nBars);
    std::vector<double> HH(nBars, std::numeric_limits<double>::lowest());
    std::vector<double> LL(nBars, std::numeric_limits<double>::max());
    std::vector<double> E(nBars), DD(nBars), PnL(nBars), trades(nBars);

    // 3) Load bars
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string date_s, time_s;
        std::string o_str, h_str, l_str, c_str, v_str;

        if (!std::getline(ss, date_s, ',') ||
            !std::getline(ss, time_s, ',') ||
            !std::getline(ss, o_str, ',') ||
            !std::getline(ss, h_str, ',') ||
            !std::getline(ss, l_str, ',') ||
            !std::getline(ss, c_str, ',') ||
            !std::getline(ss, v_str, ',')) {
            std::cerr << "Error parsing line: " << line << "\n";
            continue;
        }

        Bar bar;
        bar.t = parse_datetime(date_s, time_s);
        bar.open = std::stod(o_str);
        bar.high = std::stod(h_str);
        bar.low = std::stod(l_str);
        bar.close = std::stod(c_str);
        bar.volume = std::stod(v_str);
        bars.push_back(bar);
    }
    file.close();
    int N = static_cast<int>(bars.size());

    // 4) Determine sample indices
    std::time_t t_in_start = parse_datetime(IN_START_DATE, "00:00");
    std::time_t t_in_end = parse_datetime(IN_END_DATE, "00:00");
    std::time_t t_out_end = parse_datetime(OUT_END_DATE, "00:00");
    int in1 = 0, in2 = 0, out1 = 0, out2 = 0;
    for (int i = 0; i < N; i++) {
        if (bars[i].t >= t_in_start && in1 == 0) in1 = i;
        if (bars[i].t >= t_in_end)   in2 = i;
        if (bars[i].t >= t_in_end && out1 == 0) out1 = i;
        if (bars[i].t >= t_out_end)  out2 = i;
    }
    in1 = std::max(in1, BARS_BACK);
    in2 = std::max(in2, BARS_BACK);
    out1 = std::max(out1, BARS_BACK);
    out2 = std::max(out2, BARS_BACK);

    // 5) Rolling highs/lows using deque for O(n) performance
    if (N >= BARS_BACK && BARS_BACK >= L) {
        std::deque<int> maxDeque, minDeque;
        // Initialize window for k = BARS_BACK: window indices [BARS_BACK - L, BARS_BACK - 1]
        for (int j = BARS_BACK - L; j < BARS_BACK; j++) {
            while (!maxDeque.empty() && bars[j].high >= bars[maxDeque.back()].high)
                maxDeque.pop_back();
            maxDeque.push_back(j);

            while (!minDeque.empty() && bars[j].low <= bars[minDeque.back()].low)
                minDeque.pop_back();
            minDeque.push_back(j);
        }
        // For each k >= BARS_BACK, compute high/low for window [k-L, k)
        for (int k = BARS_BACK; k < N; k++) {
            HH[k] = bars[maxDeque.front()].high;
            LL[k] = bars[minDeque.front()].low;

            // Slide the window: remove the element leaving the window.
            int outIndex = k - L;
            if (!maxDeque.empty() && maxDeque.front() == outIndex)
                maxDeque.pop_front();
            if (!minDeque.empty() && minDeque.front() == outIndex)
                minDeque.pop_front();

            // Add the current bar into the window for the next iteration.
            if (k < N) {
                while (!maxDeque.empty() && bars[k].high >= bars[maxDeque.back()].high)
                    maxDeque.pop_back();
                maxDeque.push_back(k);
                while (!minDeque.empty() && bars[k].low <= bars[minDeque.back()].low)
                    minDeque.pop_back();
                minDeque.push_back(k);
            }
        }
    }
    else {
        // Fallback to original method (less efficient) if conditions are not met.
        for (int k = BARS_BACK; k < N; k++) {
            double hmax = std::numeric_limits<double>::lowest();
            double lmin = std::numeric_limits<double>::max();
            for (int j = k - L; j < k; j++) {
                hmax = std::max(hmax, bars[j].high);
                lmin = std::min(lmin, bars[j].low);
            }
            HH[k] = hmax;
            LL[k] = lmin;
        }
    }

    // 6) Simulate strategy
    int    position = 0;
    double equityMax = E0;
    for (int k = 0; k < N; k++) {
        E[k] = (k == 0 ? E0 : E[k - 1]);
        DD[k] = 0.0;
        PnL[k] = 0.0;
        trades[k] = 0.0;
    }

    for (int k = BARS_BACK; k < N; k++) {
        double delta = PV * (bars[k].close - bars[k - 1].close) * position;
        bool   traded = false;
        double barO = bars[k].open;
        double entryHH = HH[k];
        double entryLL = LL[k];
        double stopLong = entryHH * (1.0 - S);
        double stopShort = entryLL * (1.0 + S);

        // 6.1 Gap‐handling at OPEN
        if (position == 0) {
            if (barO >= entryHH) {
                delta = -SL_PG / 2 + PV * (bars[k].close - barO);
                position = +1;
                trades[k] = 0.5;
                traded = true;
            }
            else if (barO <= entryLL) {
                delta = -SL_PG / 2 - PV * (bars[k].close - barO);
                position = -1;
                trades[k] = 0.5;
                traded = true;
            }
        }
        else if (position == +1 && !traded) {
            if (barO <= stopLong) {
                delta = -SL_PG / 2 - PV * (barO - entryHH);
                position = 0;
                trades[k] = 0.5;
                traded = true;
            }
        }
        else if (position == -1 && !traded) {
            if (barO >= stopShort) {
                delta = -SL_PG / 2 + PV * (entryLL - barO);
                position = 0;
                trades[k] = 0.5;
                traded = true;
            }
        }

        // 6.2 Regular breakout/stop logic
        if (!traded) {
            if (position == 0) {
                bool buySignal = bars[k].high >= entryHH;
                bool sellSignal = bars[k].low <= entryLL;
                if (buySignal && !sellSignal) {
                    delta = -SL_PG / 2 + PV * (bars[k].close - entryHH);
                    position = +1;
                    trades[k] = 0.5;
                }
                else if (sellSignal && !buySignal) {
                    delta = -SL_PG / 2 - PV * (bars[k].close - entryLL);
                    position = -1;
                    trades[k] = 0.5;
                }
            }
            else if (position == +1) {
                bool revShort = bars[k].low <= entryLL;
                bool stopLOS = bars[k].low <= stopLong;
                if (revShort) {
                    delta = delta - SL_PG - 2 * PV * (bars[k].close - entryLL);
                    position = -1;
                    trades[k] = 1.0;
                }
                else if (stopLOS) {
                    delta = delta - SL_PG / 2 - PV * (bars[k].close - stopLong);
                    position = 0;
                    trades[k] = 0.5;
                }
            }
            else if (position == -1) {
                bool revLong = bars[k].high >= entryHH;
                bool stopLOS = bars[k].high >= stopShort;
                if (revLong) {
                    delta = delta - SL_PG + 2 * PV * (bars[k].close - entryHH);
                    position = +1;
                    trades[k] = 1.0;
                }
                else if (stopLOS) {
                    delta = delta - SL_PG / 2 + PV * (stopShort - bars[k].close);
                    position = 0;
                    trades[k] = 0.5;
                }
            }
        }

        // 6.3 Update equity & drawdown
        E[k] = E[k - 1] + delta;
        equityMax = std::max(equityMax, E[k]);
        DD[k] = E[k] - equityMax;
        PnL[k] = delta;
    }

    // 7) Compute & print metrics
    auto compute_stats = [&](int a, int b) {
        double profit = E[b] - E[a];
        double worstDD = 0.0, sumSq = 0.0, cnt = 0.0;
        for (int i = a; i <= b; ++i) {
            worstDD = std::min(worstDD, DD[i]);
            sumSq += PnL[i] * PnL[i];
            cnt += trades[i];
        }
        double stdev = std::sqrt(sumSq / (b - a + 1));
        return std::tuple<double, double, double, double>(profit, worstDD, stdev, cnt);
        };

    auto [pi, ddi, si, ti] = compute_stats(in1, in2);
    auto [po, ddo, so, to] = compute_stats(out1, out2);

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "In-sample:   Profit=" << pi << "  WorstDD=" << ddi << "  StDev=" << si << "  #trades=" << ti << "\n";
    std::cout << "Out-sample:  Profit=" << po << "  WorstDD=" << ddo << "  StDev=" << so << "  #trades=" << to << "\n";

    return 0;
}
