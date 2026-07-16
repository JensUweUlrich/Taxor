/*!\file SafeMap.hpp
 * \brief A mutex-protected wrapper around std::unordered_map providing coarse-grained thread-safe access.
 *
 * \details Not currently compiled into the taxor binary — this file is only #include'd by build.hpp, which is
 *          itself not #include'd by anything reachable from main.cpp as of this writing, so this class is
 *          effectively dead code too.
 */

#pragma once

#include <condition_variable>
#include <mutex>
#include <unordered_map>
#include <type_traits>

#ifndef SAFEMAP_HPP_
#define SAFEMAP_HPP_

/*!\brief Thread-safe wrapper around std::unordered_map< Key, Value >: every operation takes the internal mutex
 *        for the duration of the call, so individual operations are atomic with respect to each other, but
 *        no atomicity is provided across multiple calls (e.g. a contains() followed by a find() can race with
 *        another thread's erase() in between).
 *
 * \details \p Key and \p Value must both be move-constructible and move-assignable (enforced via the \p Enabled
 *          SFINAE parameter); this is required internally to move pairs during insertion.
 *
 * \note The two condition variables (\c cv_push, \c cv_pop) and \c push_over / \c max_size are set up as if
 *       this map were meant to act as a bounded producer/consumer queue (e.g. blocking insert() when the map is
 *       full, blocking find/erase when empty until notified), mirroring a typical "safe queue" pattern. However,
 *       in the code as written, insert()/assign() only call \c notify_one() without ever waiting on \c cv_push,
 *       and erase() only notifies \c cv_push without any consumer waiting on it — so no operation here actually
 *       blocks on the condition variables, and \c max_size / \c push_over are never enforced or read outside of
 *       set_max_size()/notify_push_over() themselves.
 */
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
        //!\brief Construct an empty map with effectively unlimited capacity (max_size = SIZE_MAX).
        SafeMap()
        : max_size( std::numeric_limits< size_t >::max() )
        {
        }


        //!\brief Construct an empty map with a given capacity limit.
        //!\param max Intended maximum size; see class-level note — this is not currently enforced by insert().
        SafeMap( size_t max )
        : max_size( max )
        {
        }

        //!\brief Update the capacity limit under the map's mutex.
        //!\param max New maximum size; see class-level note — this is not currently enforced by insert().
        void set_max_size( size_t max )
        {
            std::lock_guard< std::mutex > lock( m );
            max_size = max;
        }



        //!\brief Insert a key/value pair; a no-op if \p p.first is already present (std::unordered_map::insert
        //!       semantics — use assign() to overwrite an existing key).
        //!\param p The key/value pair to insert.
        void insert( std::pair<Key, Value> p )
        {
            // accquire mutex to modify queue
            std::unique_lock< std::mutex > lock( m );
            // insert pair in the map
            q.insert(p);
            // notify other thread that something is in the map
            cv_pop.notify_one();
        }

        //!\brief Look up (and default-construct if absent) the value for \p k, thread-safely.
        //!\param k The key to look up.
        //!\return A copy of the mapped value (default-constructed and inserted into the map if \p k was absent).
        Value operator[](const Key& k)
        {
            std::unique_lock< std::mutex > lock(m);
            return q[k];
        }

        //!\brief Insert or overwrite the value for \p p.first, thread-safely.
        //!\param p The key/value pair to insert or assign.
        void assign(std::pair<Key, Value> p)
        {
            // accquire mutex to modify queue
            std::unique_lock< std::mutex > lock(m);
            // insert pair in the map
            q.insert_or_assign(p.first, p.second);
            // notify other thread that something is in the map
            cv_pop.notify_one();
        }

        //!\brief Remove the entry for \p k, if present, thread-safely.
        //!\param k The key to remove.
        void erase( Key k)
        {
            std::unique_lock< std::mutex > lock( m );
            q.erase(k);
            cv_push.notify_one();
        }

        //!\brief Remove the entry at iterator \p t, thread-safely.
        //!\param t Iterator (obtained from a prior find()/begin()/end() call on this map) to the entry to erase.
        //!\return Iterator following the removed element, as returned by std::unordered_map::erase().
        //!\note Because \p t was obtained under a separate, already-released lock, it may be invalidated by a
        //!      concurrent modification before this call acquires the mutex; no protection against that is provided.
        typename std::unordered_map<Key, Value>::iterator erase(typename std::unordered_map<Key, Value>::iterator t)
        {
            std::unique_lock< std::mutex > lock(m);
            typename std::unordered_map<Key, Value>::iterator result = q.erase(t);
            cv_push.notify_one();
            return result;
        }

        //!\brief Mark the map as "push over" and wake all threads waiting on \c cv_pop.
        //!\note No code in this class currently waits on \c cv_pop, so this has no observable effect as written.
        void notify_push_over()
        {
            std::lock_guard< std::mutex > lock( m );
            push_over = true;
            cv_pop.notify_all();
        }

        //!\brief Current number of entries, thread-safely.
        int size()
        {
            std::lock_guard< std::mutex > lock( m );
            return q.size();
        }

        //!\brief Whether the map currently holds no entries, thread-safely.
        bool empty()
        {
            std::lock_guard< std::mutex > lock( m );
            return q.empty();
        }

        //!\brief Whether \p t is currently a key in the map, thread-safely.
        //!\param t The key to look for.
        bool contains(Key t)
        {
            std::lock_guard< std::mutex > lock(m);
            if (q.find(t) == q.end())
                return false;
            else
                return true;
        }

        //!\brief Thread-safe find(); note the returned iterator is only valid as long as no other thread
        //!       concurrently modifies the map (the lock is released before the iterator is returned to the caller).
        //!\param t The key to look for.
        //!\return Iterator to the entry, or end() if not found.
        typename std::unordered_map<Key, Value>::iterator find(Key t)
        {
            std::lock_guard< std::mutex > lock(m);
            return q.find(t);
        }

        //!\brief Thread-safe end(); see find() for the caveat about iterator validity after the lock is released.
        typename std::unordered_map<Key, Value>::iterator end()
        {
            std::lock_guard< std::mutex > lock(m);
            return q.end();
        }

        //!\brief Thread-safe begin(); see find() for the caveat about iterator validity after the lock is released.
        typename std::unordered_map<Key, Value>::iterator begin()
        {
            std::lock_guard< std::mutex > lock(m);
            return q.begin();
        }
};

#endif // SAFEMAP_HPP_
