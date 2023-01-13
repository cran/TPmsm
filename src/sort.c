
#include <stdlib.h>
#include <R_ext/Arith.h>
#include <R_ext/Boolean.h>
#include "defines.h"

static int icmp(
	int x,
	int y,
	Rboolean nalast)
{
	if (x == NA_INTEGER && y == NA_INTEGER) return 0;
	if (x == NA_INTEGER) return nalast?1:-1;
	if (y == NA_INTEGER) return nalast?-1:1;
	if (x < y)		return -1;
	if (x > y)		return 1;
	return 0;
} // icmp

static int rcmp(
	double x,
	double y,
	Rboolean nalast)
{
  int nax = R_IsNA(x) || R_IsNaN(x);
  int nay = R_IsNA(y) || R_IsNaN(y);
	if (nax && nay)	return 0;
	if (nax)		return nalast?1:-1;
	if (nay)		return nalast?-1:1;
	if (x < y)		return -1;
	if (x > y)		return 1;
	return 0;
} // rcmp

#define sub_body \
	int c = cmp(x, y, nalast); \
	if (decreasing) c = -c; \
	if (c > 0) return 1; \
	else return 0; \

static int isub(
	int x,
	int y,
	Rboolean nalast,
	Rboolean decreasing)
{
	#define cmp icmp
	sub_body
	#undef cmp
} // isub

static int rsub(
	double x,
	double y,
	Rboolean nalast,
	Rboolean decreasing)
{
	#define cmp rcmp
	sub_body
	#undef cmp
} // rsub

#define sort_body0 \
	int i, j, h; \
	for (h = 1; h <= n / 9; h = 3 * h + 1); \
	for (; h > 0; h /= 3) { \
		for (i = h; i < n; i++) { \
			v = x[i]; \
			j = i; \
			while ( j >= h && sub(x[j - h], v, nalast, decreasing) ) { \
				x[j] = x[j-h]; \
				j -= h; \
			} \
			x[j] = v; \
		} \
	} \

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Sorts vector 'x'.

Parameters:
  x[inout]          pointer to x first element.
  n[in]             length of x.
  nalast[in]        if TRUE NA values are put last.
  decreasing[in]    if TRUE sorts by descending order.

Return value:
  This function doesn't return a value.
*/

void sort_i(
	intCP x,
	int n,
	Rboolean nalast,
	Rboolean decreasing)
{
	int v;
	#define sub isub
	sort_body0
	#undef sub
} // sort_i

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Sorts vector 'x'.

Parameters:
  x[inout]          pointer to x first element.
  n[in]             length of x.
  nalast[in]        if TRUE NA values are put last.
  decreasing[in]    if TRUE sorts by descending order.

Return value:
  This function doesn't return a value.
*/

void sort_d(
	doubleCP x,
	int n,
	Rboolean nalast,
	Rboolean decreasing)
{
	double v;
	#define sub rsub
	sort_body0
	#undef sub
} // sort_d

#define sort_body \
	int i, j, h; \
	for (h = 1; h <= n / 9; h = 3 * h + 1); \
	for (; h > 0; h /= 3) { \
		for (i = h; i < n; i++) { \
			v = x[i]; \
			iv = indx[i]; \
			j = i; \
			while ( j >= h && sub(x[j - h], v, nalast, decreasing) ) { \
				x[j] = x[j-h]; \
				indx[j] = indx[j-h]; \
				j -= h; \
			} \
			x[j] = v; \
			indx[j] = iv; \
		} \
	} \

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Sorts vector 'x' with vector 'indx' alongside.

Parameters:
  x[inout]          pointer to x first element.
  indx[inout]       pointer to indx first element.
  n[in]             length of x and indx.
  nalast[in]        if TRUE NA values are put last.
  decreasing[in]    if TRUE sorts by descending order.

Return value:
  This function doesn't return a value.

Remarks:
  Vectors x and indx must have the same length.
*/

void sort_ii(
	intCP x,
	intCP indx,
	int n,
	Rboolean nalast,
	Rboolean decreasing)
{
	int v, iv;
	#define sub isub
	sort_body
	#undef sub
} // sort_ii

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Sorts vector 'x' with vector 'indx' alongside.

Parameters:
  x[inout]          pointer to x first element.
  indx[inout]       pointer to indx first element.
  n[in]             length of x and indx.
  nalast[in]        if TRUE NA values are put last.
  decreasing[in]    if TRUE sorts by descending order.

Return value:
  This function doesn't return a value.

Remarks:
  Vectors x and indx must have the same length.
*/

void sort_di(
	doubleCP x,
	intCP indx,
	int n,
	Rboolean nalast,
	Rboolean decreasing)
{
	double v;
	int iv;
	#define sub rsub
	sort_body
	#undef sub
} // sort_di

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Sorts vector 'x' with vector 'indx' alongside.

Parameters:
  x[inout]          pointer to x first element.
  indx[inout]       pointer to indx first element.
  n[in]             length of x and indx.
  nalast[in]        if TRUE NA values are put last.
  decreasing[in]    if TRUE sorts by descending order.

Return value:
  This function doesn't return a value.

Remarks:
  Vectors x and indx must have the same length.
*/

void sort_dd(
	doubleCP x,
	doubleCP indx,
	int n,
	Rboolean nalast,
	Rboolean decreasing)
{
	double v, iv;
	#define sub rsub
	sort_body
	#undef sub
} // sort_dd

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes the permutation that results by sorting time
    by increasing order.

Parameters:
  time[in]          pointer to time first element.
  index[out]        pointer to index first element.
  len[in]           length of time and index.
  nalast[in]        if TRUE NA values are put last.
  decreasing[in]    if TRUE sorts by descending order.

Return value:
  This function doesn't return a value.

Remarks:
  Vectors time and index must have the same length.
*/

void order_d(
	CdoubleCP time,
	intCP index,
	int len,
	Rboolean nalast,
	Rboolean decreasing,
	double WORK[len])
{
	register int i;
	for (i = 0; i < len; i++) WORK[i] = time[index[i]];
	sort_di(WORK, index, len, nalast, decreasing); // sort data
	return;
} // order_d

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes the permutation that results by sorting time
    by increasing order with event alongside,
    where in the subsets of constant time observations,
    event observations are sorted by decreasing order.

Parameters:
  time[in]          pointer to time first element.
  event[in]         pointer to event first element.
  index[inout]      pointer to index first element.
  len[in]           length of time, event and index.
  nalast[in]        if TRUE NA values are put last.
  decreasing0[in]   if TRUE sorts time by descending order.
  decreasing1[in]   if TRUE sorts event by descending order.

Return value:
  This function doesn't return a value.

Remarks:
  Vectors time, event and index must have the same length.
*/

void order_di(
	CdoubleCP time,
	CintCP event,
	intCP index,
	int len,
	Rboolean nalast,
	Rboolean decreasing0,
	Rboolean decreasing1,
	double WORK0[len],
	int WORK1[len])
{
	register int i;
	int j, k;
	order_d(time, index, len, nalast, decreasing0, WORK0); // compute permutation
	for (i = 0; i < len; i++) WORK1[i] = event[index[i]]; // backup vector
	i = 0; // reinitialize counter
	while (i < len) { // loop through the sample until last index is reached
		for (i++, j = 1; i < len && time[index[i]] == time[index[i-1]]; i++) { // loop through the sample until time changes or last index is reached
			j++; // count equal times
		}
		if (j > 1) { // if there are equal times
			k = i-j;
			sort_ii(&WORK1[k], &index[k], j, nalast, decreasing1); // put censored observations last
		}
	}
	return;
} // order_di

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes the permutation that results by sorting time
    by increasing order with event alongside,
    where in the subsets of constant time observations,
    event observations are sorted by decreasing order.

Parameters:
  time[in]          pointer to time first element.
  event[in]         pointer to event first element.
  index[inout]      pointer to index first element.
  len[in]           length of time, event and index.
  nalast[in]        if TRUE NA values are put last.
  decreasing0[in]   if TRUE sorts time by descending order.
  decreasing1[in]   if TRUE sorts event by descending order.

Return value:
  This function doesn't return a value.

Remarks:
  Vectors time, event and index must have the same length.
*/

void order_dd(
	CdoubleCP time,
	CdoubleCP event,
	intCP index,
	int len,
	Rboolean nalast,
	Rboolean decreasing0,
	Rboolean decreasing1,
	double WORK0[len],
	double WORK1[len])
{
	register int i;
	int j, k;
	order_d(time, index, len, nalast, decreasing0, WORK0); // compute permutation
	for (i = 0; i < len; i++) WORK1[i] = event[index[i]]; // backup vector
	i = 0; // reinitialize counter
	while (i < len) { // loop through the sample until last index is reached
		for (i++, j = 1; i < len && time[index[i]] == time[index[i-1]]; i++) { // loop through the sample until time changes or last index is reached
			j++; // count equal times
		}
		if (j > 1) { // if there are equal times
			k = i-j;
			sort_di(&WORK1[k], &index[k], j, nalast, decreasing1); // put censored observations last
		}
	}
	return;
} // order_dd
