#pragma once

// LRU cache primitives used by CASContext.
// Extracted from symbolic.hpp to keep that file within the 500-line anti-monolith limit.

#include <cstddef>
#include <cstdint>
#include <list>
#include <optional>
#include <unordered_map>
#include <utility>

namespace cas::symbolic {

struct CacheMetrics {
    std::uint64_t hits{0};
    std::uint64_t misses{0};
    std::uint64_t evictions{0};
};

template <typename Key, typename Value, typename Hash = std::hash<Key>, typename Equal = std::equal_to<Key>>
class CacheContainer {
public:
    using ListType = std::list<Key>;
    using MapType = std::unordered_map<Key, std::pair<Value, typename ListType::iterator>, Hash, Equal>;

    explicit CacheContainer(std::size_t max_size = 1000) : max_size_(max_size) {}

    void set_max_size(std::size_t size) {
        max_size_ = size;
        evict_if_needed();
    }

    [[nodiscard]] std::size_t max_size() const noexcept { return max_size_; }
    [[nodiscard]] std::size_t size() const noexcept { return map_.size(); }

    [[nodiscard]] std::optional<Value> get(const Key& key) {
        auto it = map_.find(key);
        if (it == map_.end()) {
            metrics_.misses++;
            return std::nullopt;
        }
        metrics_.hits++;
        list_.splice(list_.begin(), list_, it->second.second);
        return it->second.first;
    }

    void put(const Key& key, Value value) {
        if (max_size_ == 0) return;

        auto it = map_.find(key);
        if (it != map_.end()) {
            it->second.first = value;
            list_.splice(list_.begin(), list_, it->second.second);
            return;
        }

        list_.push_front(key);
        map_[key] = {value, list_.begin()};
        evict_if_needed();
    }

    void clear() noexcept {
        map_.clear();
        list_.clear();
    }

    [[nodiscard]] CacheMetrics& metrics() noexcept { return metrics_; }
    [[nodiscard]] const CacheMetrics& metrics() const noexcept { return metrics_; }
    void reset_metrics() noexcept { metrics_ = {}; }

    [[nodiscard]] auto begin() { return map_.begin(); }
    [[nodiscard]] auto end() { return map_.end(); }
    [[nodiscard]] auto begin() const { return map_.begin(); }
    [[nodiscard]] auto end() const { return map_.end(); }
    [[nodiscard]] bool empty() const noexcept { return map_.empty(); }

private:
    void evict_if_needed() {
        while (map_.size() > max_size_ && !list_.empty()) {
            Key last = list_.back();
            list_.pop_back();
            map_.erase(last);
            metrics_.evictions++;
        }
    }

    std::size_t max_size_;
    MapType map_;
    ListType list_;
    CacheMetrics metrics_;
};

}  // namespace cas::symbolic
