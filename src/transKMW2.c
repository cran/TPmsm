/*************************************/
/*** KMW2 TRANSITION PROBABILITIES ***/
/*************************************/

#include "defines.h"
#include "get.h"

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Computes the transition probabilities:
		p11(s,t) = P(Z>t|Z>s) = {1-P(Z<=t)}/{1-P(Z<=s)}
		p12(s,t) = 1-p11(s,t)-p13(s,t)
		p13(s,t) = P(T<=t|Z>s) = P(Z>s,T<=t)/{1-P(Z<=s)}
		p23(s,t) = P(T<=t|Z<=s,T>s) = P(Z<=s,s<T<=t)/{P(Z<=s)-P(T<=s)}

Parameters:
	len[in]			pointer to length of T1, E1, S and E.
	T1[in]			pointer to T1 first element.
	E1[in]			pointer to E1 first element.
	S[in]			pointer to S first element.
	E[in]			pointer to E first element.
	index0[in]		pointer to index0 first element.
	index1[in]		pointer to index1 first element.
	nt[in]			pointer to length of UT and number of rows of P.
	UT[in]			pointer to unique times vector.
	nb[in]			pointer to number of rows of P.
	P[out]			pointer to a (nb)x(nt)x4 probability array.
	b[in]			pointer to row index.

Return value:
	This function doesn't return a value.

Remarks:
	Vector index0 must indicate the permutation of vector T1
		sorted by ascending order.
	Vector index1 must indicate the permutation of vector S
		sorted by ascending order.
	Vectors T1, E1, S and E must have the same length.
*/

void transKMW3I(
	CintCP len,
	Cdouble T1[*len],
	Cint E1[*len],
	Cdouble S[*len],
	Cint E[*len],
	Cint index0[*len],
	Cint index1[*len],
	CintCP nt,
	Cdouble UT[*nt],
	CintCP nb,
	double P[*nb*(*nt)*4],
	CintCP b)
{
	register int i;
	int j = 0, e;
	double aux[4], p[2];
	getIndexI(T1, index0, &UT[0], len, &j, &e); // determine first index
	for (aux[0] = 1, aux[3] = 1; j < e; j++) { // loop through the sample until last index is reached
		aux[2] = (double)E1[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] -= aux[2]; // compute probability
	}
	getIndexI(T1, index0, &UT[*nt-1], len, &j, &e); // determine last index
	for (i = 0; j < e; j++) { // loop through the sample until last index is reached
		while (T1[index0[j]] > UT[i]) P[*b+*nb*i++] = aux[3];
		aux[2] = (double)E1[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] -= aux[2]; // compute probability
	}
	for (; i < *nt; i++) P[*b+*nb*i] = aux[3];
	j = 0;
	getIndexI(S, index1, &UT[0], len, &j, &e); // determine first index
	for (aux[0] = 1, aux[3] = 0; j < e; j++) { // loop through the sample until last index is reached
		aux[2] = (double)E[index1[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] += aux[2]; // compute probability
	}
	getIndexI(S, index1, &UT[*nt-1], len, &j, &e); // determine last index
	for (i = 0, p[0] = 0, p[1] = 0; j < e; j++) { // loop through the sample until last index is reached
		while (S[index1[j]] > UT[i]) {
			P[*b+*nb*(i+*nt*2)] = p[0];
			P[*b+*nb*(i+*nt*3)] = p[1];
			i++;
		}
		aux[2] = (double)E[index1[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		if (T1[index1[j]] <= UT[0]) p[1] += aux[2]; // compute probability
		else p[0] += aux[2]; // compute probability
	}
	for (; i < *nt; i++) {
		P[*b+*nb*(i+*nt*2)] = p[0];
		P[*b+*nb*(i+*nt*3)] = p[1];
	}
	for (i = *nt-1; i >= 0; i--) {
		P[*b+*nb*(i+*nt*2)] /= P[*b]; // compute and save p13(s,t)
		P[*b+*nb*(i+*nt*3)] /= 1-P[*b]-aux[3]; // compute and save p23(s,t)
		if (P[*b+*nb*(i+*nt*3)] > 1) P[*b+*nb*(i+*nt*3)] = 1;
		P[*b+*nb*i] /= P[*b]; // compute and save p11(s,t)
		P[*b+*nb*(i+*nt)] = 1-P[*b+*nb*i]-P[*b+*nb*(i+*nt*2)]; // compute and save p12(s,t)
		if (P[*b+*nb*(i+*nt)] < 0) {
			P[*b+*nb*(i+*nt*2)] = 1-P[*b+*nb*i]; // compute and save p13(s,t)
			P[*b+*nb*(i+*nt)] = 0; // compute and save p12(s,t)
		}
	}
	return;
} // transKMW3I

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Computes the transition probabilities:
		p11(s,t) = P(Z>t|Z>s) = P(Z>t)/P(Z>s)
		p12(s,t) = 1-p11(s,t)-p13(s,t)
		p13(s,t) = P(T<=t|Z>s) = P(Z>s,T<=t)/P(Z>s)
		p23(s,t) = P(T<=t|Z<=s,T>s) = P(Z<=s,s<T<=t)/P(Z<=s,T>s)

Parameters:
	len[in]			pointer to length of T1, E1, S and E.
	T1[in]			pointer to T1 first element.
	E1[in]			pointer to E1 first element.
	S[in]			pointer to S first element.
	E[in]			pointer to E first element.
	index0[in]		pointer to index0 first element.
	index1[in]		pointer to index1 first element.
	nt[in]			pointer to length of UT and number of rows of P.
	UT[in]			pointer to unique times vector.
	nb[in]			pointer to number of rows of P.
	P[out]			pointer to a (nb)x(nt)x4 probability array.
	b[in]			pointer to row index.

Return value:
	This function doesn't return a value.

Remarks:
	Vector index0 must indicate the permutation of vector T1
		sorted by ascending order.
	Vector index1 must indicate the permutation of vector S
		sorted by ascending order.
	Vectors T1, E1, S and E must have the same length.
*/

void transKMW4I(
	CintCP len,
	Cdouble T1[*len],
	Cint E1[*len],
	Cdouble S[*len],
	Cint E[*len],
	Cint index0[*len],
	Cint index1[*len],
	CintCP nt,
	Cdouble UT[*nt],
	CintCP nb,
	double P[*nb*(*nt)*4],
	CintCP b)
{
	register int i;
	int j = 0, e;
	double aux[4], p[2];
	getIndexI(T1, index0, &UT[0], len, &j, &e); // determine first index
	for (aux[0] = 1, aux[3] = 0; j < e; j++) { // loop through the sample until last index is reached
		aux[2] = (double)E1[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] += aux[2]; // compute probability
	}
	getIndexI(T1, index0, &UT[*nt-1], len, &j, &e); // determine last index
	for (i = 0; j < e; j++) { // loop through the sample until last index is reached
		while (T1[index0[j]] > UT[i]) P[*b+*nb*i++] = aux[3];
		aux[2] = (double)E1[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] += aux[2]; // compute probability
	}
	for (; i < *nt; i++) P[*b+*nb*i] = aux[3];
	for (; j < *len; j++) { // loop through the sample until last index is reached
		aux[2] = (double)E1[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] += aux[2]; // compute probability
	}
	for (i = 0; i < *nt; i++) P[*b+*nb*i] = aux[3]-P[*b+*nb*i];
	j = 0;
	getIndexI(S, index1, &UT[0], len, &j, &e); // determine first index
	for (aux[0] = 1; j < e; j++) { // loop through the sample until last index is reached
		aux[1] = 1-(double)E[index1[j]]/(*len-j); // compute needed factor
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
	}
	getIndexI(S, index1, &UT[*nt-1], len, &j, &e); // determine last index
	for (i = 0, p[0] = 0, p[1] = 0; j < e; j++) { // loop through the sample until last index is reached
		while (S[index1[j]] > UT[i]) {
			P[*b+*nb*(i+*nt*2)] = p[0]; // save probability
			P[*b+*nb*(i+*nt*3)] = p[1]; // save probability
			i++;
		}
		aux[2] = (double)E[index1[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		if (T1[index1[j]] <= UT[0]) p[1] += aux[2]; // compute probability
		else p[0] += aux[2]; // compute probability
	}
	for (; i < *nt; i++) {
		P[*b+*nb*(i+*nt*2)] = p[0];
		P[*b+*nb*(i+*nt*3)] = p[1];
	}
	for (; j < *len; j++) { // loop through the sample until last index is reached
		aux[2] = (double)E[index1[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		if (T1[index1[j]] <= UT[0]) p[1] += aux[2]; // compute probability
	}
	for (i = *nt-1; i >= 0; i--) {
		P[*b+*nb*(i+*nt*2)] /= P[*b]; // compute and save p13(s,t)
		P[*b+*nb*i] /= P[*b]; // compute and save p11(s,t)
		P[*b+*nb*(i+*nt)] = 1-P[*b+*nb*i]-P[*b+*nb*(i+*nt*2)]; // compute and save p12(s,t)
		if (P[*b+*nb*(i+*nt)] < 0) {
			P[*b+*nb*(i+*nt*2)] = 1-P[*b+*nb*i]; // compute and save p13(s,t)
			P[*b+*nb*(i+*nt)] = 0; // compute and save p12(s,t)
		}
		P[*b+*nb*(i+*nt*3)] /= p[1]; // compute and save p23(s,t)
	}
	return;
} // transKMW4I
