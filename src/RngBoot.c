
#include "defines.h"
#include "RngStream.h"

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes bootstrap index.

Parameters:
  g[in]             pointer to RngStream_InfoState.
  len[in]           pointer to length of index.
  index[out]        pointer to index first element.

Return value:
  This function doesn't return a value.
*/

void boot_i(
	RngStream g,
	CintCP len,
	int index[*len])
{
	register int i;
	for (i = 0; i < *len; i++) {
		index[i] = (int)RngStream_RandInt(g, 0, (long)*len-1);
	}
	return;
} // boot_i

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes bootstrap indexes.

Parameters:
  g[in]             pointer to RngStream_InfoState.
  len[in]           pointer to length of index0 and index1.
  index0[out]       pointer to index0 first element.
  index1[out]       pointer to index1 first element.

Return value:
  This function doesn't return a value.

Remarks:
  Vectors index0 and index1 must have the same length.
*/

void boot_ii(
	RngStream g,
	CintCP len,
	int index0[*len],
	int index1[*len])
{
	register int i;
	for (i = 0; i < *len; i++) {
		index0[i] = (int)RngStream_RandInt(g, 0, (long)*len-1);
		index1[i] = index0[i];
	}
	return;
} // boot_ii

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Initializes index without bootstrap.

Parameters:
  len[in]           pointer to length of index.
  index[out]        pointer to index first element.

Return value:
  This function doesn't return a value.
*/

void indx_i(
  CintCP len,
  int index[*len])
{
	register int i;
	for (i = 0; i < *len; i++) index[i] = i;
	return;
} // indx_i

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Initializes indexes without bootstrap.

Parameters:
  len[in]           pointer to length of index0 and index1.
  index0[out]       pointer to index0 first element.
  index1[out]       pointer to index1 first element.

Return value:
  This function doesn't return a value.

Remarks:
  Vectors index0 and index1 must have the same length.
*/

void indx_ii(
  CintCP len,
  int index0[*len],
  int index1[*len])
{
	register int i;
	for (i = 0; i < *len; i++) {
		index0[i] = i;
		index1[i] = i;
	}
	return;
} // indx_ii
