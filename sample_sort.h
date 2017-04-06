template<typename T, typename Compare>
void parallel_sample_sort( T* xs, T* xe, Compare cmp) {
	if( xe-xs<=SAMPLE_SORT_CUT_OFF ) {
		//parallel_quicksort(xs,xe);
		std::sort(xs, xe, cmp);
	} else {
		size_t m = choose_number_of_bins(xe-xs);
		size_t tally[M_MAX][M_MAX];
		T* y = new T[xe-xs];
		bin(xs, xe, m, y, tally, cmp);
		repack_and_subsort(xs, xe, m, y, tally, cmp);
		delete[] y;
	}
}
template<typename T, typename Compare>
bool parallel_sample_sort_sampling( T* xs, T* xe, Compare cmp ) {
	size_t m = choose_number_of_bins(xe-xs);
	size_t tally[M_MAX][M_MAX];
	T* y = new T[xe-xs];
	bin(xs, xe, m, y, tally, cmp);
	long max=0; long min=0;
	for(int i = 0; i< m; ++i)
	{
		min = tally[i][m-1] < min? tally[i][m-1] : min;
		max = tally[i][m-1] > max? tally[i][m-1] : max;
	}
	delete[] y;
	return (max-min)>min;
}
