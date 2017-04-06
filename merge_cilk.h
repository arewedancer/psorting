// merge sequences [xs,xe) and [ys,ye) to output [zs,(xe-xs)+(ye-ys)
template<typename T, typename Compare>
void parallel_merge_cilk( T* xs, T* xe, T* ys, T* ye, T* zs, Compare cmp ) {
    const size_t MERGE_CUT_OFF = 2000;
    if( xe-xs + ye-ys <= MERGE_CUT_OFF ) {
        serial_merge(xs,xe,ys,ye,zs,cmp);
    } else {
        T *xm, *ym;
        if( xe-xs < ye-ys  ) {
            ym = ys+(ye-ys)/2;
            xm = std::upper_bound(xs,xe,*ym,cmp);
        } else {
            xm = xs+(xe-xs)/2;
            ym = std::lower_bound(ys,ye,*xm,cmp);
        }
        T* zm = zs + (xm-xs) + (ym-ys);
        cilk_spawn parallel_merge_cilk( xs, xm, ys, ym, zs, cmp );
        /*nospawn*/parallel_merge_cilk( xm, xe, ym, ye, zm, cmp );
        // implicit cilk_sync
    }
}
