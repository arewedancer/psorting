#include <cilk/cilk_api.h> 
// sorts [xs,xe).  zs[0:xe-xs) is temporary buffer supplied by caller.
// result is in [xs,xe) if inplace==true, otherwise in zs[0:xe-xs)
template<typename T, typename Compare>
void parallel_merge_sort_cilk( T* xs, T* xe, T* zs, bool inplace, Compare cmp ) {
	const size_t SORT_CUT_OFF = 500;
	if( xe-xs<=SORT_CUT_OFF ) {
		std::sort( xs, xe, cmp );
		if( !inplace ) 
			std::move( xs, xe, zs );
	} else {
		T* xm = xs + (xe-xs)/2;
		T* zm = zs + (xm-xs);
		T* ze = zs + (xe-xs);
		cilk_spawn parallel_merge_sort_cilk( xs, xm, zs, !inplace, cmp );
		/*nospawn*/parallel_merge_sort_cilk( xm, xe, zm, !inplace, cmp );
		cilk_sync;
		if( inplace )
			parallel_merge_cilk( zs, zm, zm, ze, xs, cmp );
		else
			parallel_merge_cilk( xs, xm, xm, xe, zs, cmp );
	}
}

template<typename T, typename Compare>
void parallel_merge_sort_cilkplus( T* xs, T* xe, Compare cmp, int max_thread=0 ) { 
    T* zs = new T[xe-xs];
		if (max_thread != 0)
			__cilkrts_set_param("nworkers", max_thread);

    parallel_merge_sort_cilk( xs, xe, zs, true, cmp );
    delete[] zs; 
}

