
void sort_i(intCP x, int n, Rboolean nalast, Rboolean decreasing);
void sort_d(doubleCP x, int n, Rboolean nalast, Rboolean decreasing);
void sort_ii(intCP x, intCP indx, int n, Rboolean nalast, Rboolean decreasing);
void sort_di(doubleCP x, intCP indx, int n, Rboolean nalast, Rboolean decreasing);
void sort_dd(doubleCP x, doubleCP indx, int n, Rboolean nalast, Rboolean decreasing);
void order_d(CdoubleCP time, intCP index, int len, Rboolean nalast, Rboolean decreasing, double WORK[len]);
void order_di(CdoubleCP time, CintCP event, intCP index, int len, Rboolean nalast, Rboolean decreasing0, Rboolean decreasing1, double WORK0[len], int WORK1[len]);
void order_dd(CdoubleCP time, CdoubleCP event, intCP index, int len, Rboolean nalast, Rboolean decreasing0, Rboolean decreasing1, double WORK0[len], double WORK1[len]);
