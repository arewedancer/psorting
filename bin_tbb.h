template<typename T, typename Compare>
void bin( T* xs, T* xe, size_t m, T* y, size_t tally[M_MAX][M_MAX], Compare cmp) {
    T tree[M_MAX-1];
    build_sample_tree( xs, xe, tree, m, cmp);
		//exit(1);

    size_t block_size = ((xe-xs)+m-1)/m;
    bindex_type* bindex = new bindex_type[xe-xs];
    tbb::parallel_for( size_t(0), m, [=,&tree](size_t i) { 
        size_t js = i*block_size;
        size_t je = std::min( js+block_size, size_t(xe-xs) );

        // Map keys to bins
        size_t freq[M_MAX];
        map_keys_to_bins( xs+js, je-js, tree, m, bindex+js, freq, cmp );

        // Compute where each bin starts
        T* dst[M_MAX];
        size_t s = 0;
        for( size_t j=0; j<m; ++j ) {
            dst[j] = y+js+s;
            s += freq[j];
            tally[i][j] = s;
        }

        // Scatter keys into their respective bins
        for( size_t j=js; j<je; ++j )
            *dst[bindex[j]]++ = std::move(xs[j]);
    });
    delete[] bindex;
}
