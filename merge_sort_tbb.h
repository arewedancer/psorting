#include "tbb/task_scheduler_init.h"
// sorts [xs,xe).  zs[0:xe-xs) is temporary buffer supplied by caller.
// result is in [xs,xe) if inplace==true, otherwise in zs[0:xe-xs)
template<typename T, typename Compare>
void parallel_merge_sort( T* xs, T* xe, T* zs, bool inplace, Compare cmp ) {
    const size_t SORT_CUT_OFF = 500;
    if( xe-xs<=SORT_CUT_OFF ) {
        std::sort( xs, xe, cmp );
        if( !inplace ) 
				{
            std::move( xs, xe, zs );
				};
    } else {
       T* xm = xs + (xe-xs)/2;
       T* zm = zs + (xm-xs);
       T* ze = zs + (xe-xs);
       tbb::parallel_invoke( [=]{parallel_merge_sort( xs, xm, zs, !inplace, cmp );},
                             [=]{parallel_merge_sort( xm, xe, zm, !inplace, cmp );} );
       if( inplace )
			 {
           parallel_merge( zs, zm, zm, ze, xs, cmp );
			 }
       else
			 {
           parallel_merge( xs, xm, xm, xe, zs, cmp );
			 }
   }
}
template<typename T, typename Compare>
void parallel_merge_sort_tbb( T* xs, T* xe, Compare cmp, int max_thread=0 ) {
    T* zs = new T[xe-xs];
		if (max_thread != 0)
			 tbb::task_scheduler_init init(max_thread);
    parallel_merge_sort( xs, xe, zs, true, cmp );
    delete[] zs;
}
