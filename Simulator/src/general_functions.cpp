#include "general_functions.h"

//sort fragment
void merge_general(void** sortlist, double* sortkey, long p, long q, long r, double* mergeSort_Larray, double* mergeSort_Rarray, void** mergeSort_LorderedList, void** mergeSort_RorderedList)
{
	long n1, n2, i, j, k;

	n1 = q - p + 1;
	n2 = r - q;

	for (i = 1; i <= n1; i++)
	{
		mergeSort_Larray[i] = sortkey[p + i - 1];
		mergeSort_LorderedList[i] = sortlist[p + i - 1];
	}
	for (j = 1; j <= n2; j++)
	{
		mergeSort_Rarray[j] = sortkey[q + j];
		mergeSort_RorderedList[j] = sortlist[q + j];
	}

	mergeSort_Larray[n1 + 1] = MAX_NUMBER;
	mergeSort_Rarray[n2 + 1] = MAX_NUMBER;

	i = 1;
	j = 1;

	for (k = p; k <= r; k++)
	{
		if (mergeSort_Larray[i] <= mergeSort_Rarray[j])
		{
			sortkey[k] = mergeSort_Larray[i];
			sortlist[k] = mergeSort_LorderedList[i];

			i++;
		} 
		else
		{
			sortkey[k] = mergeSort_Rarray[j];
			sortlist[k] = mergeSort_RorderedList[j];

			j++;
		}
	}

	return;
}


void mergeSort_general(void** sortlist, double* sortkey, long sortList_size)
{
	if (sortList_size <= 0)
		return;

	//non-recursive merge sort for sorting junctions
	long m, n, i, r;
	m = 1;
	n = sortList_size;

	double* mergeSort_Larray = new double [sortList_size + 2];
	double* mergeSort_Rarray = new double [sortList_size + 2];
	void** mergeSort_LorderedList = new void* [sortList_size + 2];
	void** mergeSort_RorderedList = new void* [sortList_size + 2];

	while (m <= n)
	{
		i = 1;
		while (i <= n - m)
		{
			r = (i + 2 * m - 1) < n ? (i + 2 * m - 1) : n;
			merge_general(sortlist, sortkey, i, i + m - 1, r, mergeSort_Larray, mergeSort_Rarray, mergeSort_LorderedList, mergeSort_RorderedList);
			i = i + 2 * m;
		}

		m = m * 2;
	}

	delete [] mergeSort_Larray;
	delete [] mergeSort_Rarray;
	delete [] mergeSort_LorderedList;
	delete [] mergeSort_RorderedList;

	return;
}

