
#include <stdlib.h>
#include "defines.h"

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Computes bootstrap index.

Parameters:
	index[out]		pointer to index first element.
	len[in]			pointer to length of index.

Return value:
	This function doesn't return a value.
*/

void boot_i(
	intCP index,
	CintCP len)
{
	register int i;
	for (i = 0; i < *len; i++) index[i] = *len * (double)rand()/(RAND_MAX+1.0);
	return;
} // boot_i

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Computes bootstrap indexes.

Parameters:
	index0[out]		pointer to index0 first element.
	index1[out]		pointer to index1 first element.
	len[in]			pointer to length of index0 and index1.

Return value:
	This function doesn't return a value.

Remarks:
	Vectors index0 and index1 must have the same length.
*/

void boot_ii(
	intCP index0,
	intCP index1,
	CintCP len)
{
	register int i;
	for (i = 0; i < *len; i++) {
		index0[i] = *len * (double)rand()/(RAND_MAX+1.0);
		index1[i] = index0[i];
	}
	return;
} // boot_ii

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Initializes index without bootstrap.

Parameters:
	index[out]		pointer to index first element.
	len[in]			pointer to length of index.

Return value:
	This function doesn't return a value.
*/

void indx_i(
	intCP index,
	CintCP len)
{
	register int i;
	for (i = 0; i < *len; i++) index[i] = i;
	return;
} // indx_i

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Initializes indexes without bootstrap.

Parameters:
	index0[out]		pointer to index0 first element.
	index1[out]		pointer to index1 first element.
	len[in]			pointer to length of index0 and index1.

Return value:
	This function doesn't return a value.

Remarks:
	Vectors index0 and index1 must have the same length.
*/

void indx_ii(
	intCP index0,
	intCP index1,
	CintCP len)
{
	register int i;
	for (i = 0; i < *len; i++) {
		index0[i] = i;
		index1[i] = i;
	}
	return;
} // indx_ii
