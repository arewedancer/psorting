template<typename T, typename Compare>
void serial_merge( T* xs, T* xe, T* ys, T* ye, T* zs, Compare cmp) {
    while( xs!=xe && ys!=ye ) {
        //bool which = ascending ? *ys < *xs : *ys > *xs;
        bool which = cmp(*ys, *xs);
        *zs++ = std::move(which ? *ys++ : *xs++);
    }
    std::move( xs, xe, zs );
    std::move( ys, ye, zs );
}
