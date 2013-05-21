/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/blitz.hpp>
#include <libmpdata++/concurr/detail/sharedmem.hpp>
#include <libmpdata++/solvers/detail/monitor.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      constexpr int max(const int a, const int b)
      {
        return a > b ? a : b;
      }

      template <typename real_t_, int n_dims_, int n_eqs_, int n_tlev_, int halo_>
      class solver_common
      {
	public:

        // using enums as "static const int" would need instantiation
        enum { halo = halo_ }; 
        enum { n_dims = n_dims_ };
        enum { n_eqs = n_eqs_ };
        enum { n_tlev = n_tlev_ };

        typedef real_t_ real_t;
        typedef blitz::Array<real_t_, n_dims_> arr_t;
        typedef arrvec_t<arr_t> arrvec_t;

	void cycle_all()
	{ 
	  for (int e = 0; e < n_eqs; ++e) cycle(e);
          this->mem->cycle();
	}

	protected: 

        std::vector<int> n;

        typedef concurr::detail::sharedmem<real_t, n_dims, n_eqs, n_tlev> mem_t;
	mem_t *mem;

	// helper methods invoked by solve()
	virtual void advop(int e) = 0;
	void advop_all()
	{
	  for (int e = 0; e < n_eqs; ++e) advop(e);
	}

	void cycle(int e) 
	{ 
	  n[e] = (n[e] + 1) % n_tlev - n_tlev;  // TODO: - n_tlev not needed?
	}

	virtual void xchng(int e, int l = 0) = 0; // TODO: make l -> -l
	void xchng_all() 
	{   
	  for (int e = 0; e < n_eqs; ++e) xchng(e);
	}

        // hook methods to be overrided in inheriting objects
        // the overriding methods are expected to call parent_t::hook_...()

        private:
      
        bool 
          hook_ante_step_called = false, 
          hook_post_step_called = false, 
          hook_ante_loop_called = false, 
          hook_post_loop_called = false;

	public:

        virtual void hook_ante_step() 
        { 
#if !defined(NDEBUG)
          hook_ante_step_called = true;
#endif
        }

        virtual void hook_post_step() 
        {
#if !defined(NDEBUG)
          hook_post_step_called = true;
#endif
        }

        virtual void hook_ante_loop(const int nt) 
        {
#if !defined(NDEBUG)
          hook_ante_loop_called = true;
#endif
        }

        virtual void hook_post_loop() 
        {
#if !defined(NDEBUG)
          hook_post_loop_called = true;
#endif
        }

	// ctor
	solver_common(mem_t *mem) :
	  n(n_eqs, 0), 
          mem(mem)
	{ }

        // dtor
        virtual ~solver_common()
        {
#if !defined(NDEBUG)
          assert(hook_ante_step_called && "any overriding hook_ante_step() must call parent_t::hook_ante_step()");
          assert(hook_post_step_called && "any overriding hook_post_step() must call parent_t::hook_post_step()");
          assert(hook_ante_loop_called && "any overriding hook_ante_loop() must call parent_t::hook_ante_loop()");
          assert(hook_post_loop_called && "any overriding hook_post_loop() must call parent_t::hook_post_loop()");
#endif
        }

	virtual void solve(const int nt) final
	{   
          hook_ante_loop(nt);
	  for (int t = 0; t < nt; ++t) 
	  {   
            monitor(float(t) / nt);
            hook_ante_step();
	    xchng_all();
	    advop_all();
	    cycle_all();
            hook_post_step();
	  }   
          hook_post_loop();
        }

	// psi getter
	arr_t state(int e, int add = 0)
	{
	  return this->mem->psi[e][this->n[e] + add];
	}
      };
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx