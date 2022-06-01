/*
  Project Name	: KmerGO
  Version		: 2.0.0
  Author		: Qi Chen
  Date			: 2020-08-03
*/

#include "KmerGO.h"

/* Constructor */
KmerLoserTree::KmerLoserTree(uint32 c, uint32 k)
{
	capacity = c;
	CKmerAPI kmer_obj(k);
	for (uint32 i = 0; i < c; i++)
	{
		kmer_list.push_back(kmer_obj);
		counter_list.push_back(0);
	}
	kmer_list.push_back(kmer_obj);
}

/* Destructor */
KmerLoserTree::~KmerLoserTree()
{
	;
}

/* Swap two variables */
void KmerLoserTree::swap(uint32* a, uint32* b)
{
	int tmp;
	tmp = *a;
	*a = *b;
	*b = tmp;
}

/* Adjust the loser tree and s is the index */
void KmerLoserTree::adjust(uint32 s)
{
	int t = (s + capacity) / 2;		//get the father node of s
	while (t > 0) {					//compare with the parent node until the root node of the loser tree
		if (!kmer_list[s].operator<(kmer_list[ls[t]]) && !kmer_list[s].operator==(kmer_list[ls[t]]))
		{
			swap(&s, &ls[t]);		//swap s and ls[t]
		}
		t /= 2;
	}
	ls[0] = s;						//the last winner is in ls[0]
}

/* Build a loser tree */
void KmerLoserTree::build()
{
	for (uint32 i = 0; i < capacity; i++) ls.push_back(capacity);
	for (int i = capacity - 1; i >= 0; i--) adjust(i);	//adjust the position
}

/* Get the minimum position*/
uint32 KmerLoserTree::getMIN()
{
	return ls[0];
}
