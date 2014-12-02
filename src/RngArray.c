
#include <stddef.h>
#include <Rdefines.h>
#include "defines.h"
#include "RngStream.h"

RngStream *RngArray = NULL;

void RngArray_CreateStream(
	CintCP len,
	RngStream RngArray[*len])
{
	register int i;
	for (i = 0; i < *len; i++) {
		RngArray[i] = RngStream_CreateStream("");
	}
	return;
} // RngArray_CreateStream

void RngArray_DeleteStream(
	CintCP len,
	RngStream RngArray[*len])
{
	register int i;
	for (i = 0; i < *len; i++) {
		RngStream_DeleteStream(&RngArray[i]);
	}
	return;
} // RngArray_DeleteStream
