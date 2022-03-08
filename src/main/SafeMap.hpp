#pragma once

#include <condition_variable>
#include <mutex>
#include <unordered_map>
#include <type_traits>

#ifndef SAFEMAP_HPP_
#define SAFEMAP_HPP_

template < class Key, class Value,
           class Enabled = std::enable_if_t< std::is_move_assignable< Key >::value && //
                                             std::is_move_constructible< Key >::value &&
                                             std::is_move_assignable< Value >::value && 
                                             std::is_move_constructible< Value >::value> >
class SafeMap
{

    private:
        std::unordered_map< Key , Value > q;
        std::mutex      m;
        // Very high max_size - default to std::numeric_limits< size_t >::max() - virtually no limit on size of the map
        size_t                  max_size;
        bool                    push_over = false;
        std::condition_variable cv_push;
        std::condition_variable cv_pop;

    public:
        SafeMap()
        : max_size( std::numeric_limits< size_t >::max() )
        {
        }


        SafeMap( size_t max )
        : max_size( max )
        {
        }

        void set_max_size( size_t max )
        {
            std::lock_guard< std::mutex > lock( m );
            max_size = max;
        }

    

        void insert( std::pair<Key, Value> p )
        {
            // accquire mutex to modify queue
            std::unique_lock< std::mutex > lock( m );
            // insert pair in the map
            q.insert(p);
            // notify other thread that something is in the map
            cv_pop.notify_one();
        }

        Value operator[](const Key& k)
        {
            std::unique_lock< std::mutex > lock(m);
            return q[k];
        }

        void assign(std::pair<Key, Value> p)
        {
            // accquire mutex to modify queue
            std::unique_lock< std::mutex > lock(m);
            // insert pair in the map
            q.insert_or_assign(p.first, p.second);
            // notify other thread that something is in the map
            cv_pop.notify_one();
        }

        void erase( Key k)
        {
            std::unique_lock< std::mutex > lock( m );
            q.erase(k);
            cv_push.notify_one();
        }

        typename std::unordered_map<Key, Value>::iterator erase(typename std::unordered_map<Key, Value>::iterator t)
        {
            std::unique_lock< std::mutex > lock(m);
            typename std::unordered_map<Key, Value>::iterator result = q.erase(t);
            cv_push.notify_one();
            return result;
        }

        void notify_push_over()
        {
            std::lock_guard< std::mutex > lock( m );
            push_over = true;
            cv_pop.notify_all();
        }

        int size()
        {
            std::lock_guard< std::mutex > lock( m );
            return q.size();
        }

        bool empty()
        {
            std::lock_guard< std::mutex > lock( m );
            return q.empty();
        }

        bool contains(Key t)
        {
            std::lock_guard< std::mutex > lock(m);
            if (q.find(t) == q.end())
                return false;
            else
                return true;
        }

        typename std::unordered_map<Key, Value>::iterator find(Key t)
        {
            std::lock_guard< std::mutex > lock(m);
            return q.find(t);
        }

        typename std::unordered_map<Key, Value>::iterator end()
        {
            std::lock_guard< std::mutex > lock(m);
            return q.end();
        }

        typename std::unordered_map<Key, Value>::iterator begin()
        {
            std::lock_guard< std::mutex > lock(m);
            return q.begin();
        }
};

#endif // SAFEMAP_HPP_
