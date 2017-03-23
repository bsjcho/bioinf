package mdp

import (
	"math"

	bio "github.com/bsjcho/bioinf"
	"github.com/bsjcho/nd"
)

type multiDP struct {
	seqs   []*bio.Sequence // list of sequences
	table  *nd.Array       // dp table to store optimal scores
	cached *nd.Array       // to determine if an optimal score has already been
	// calculated. necessary for memoization since scores can be 0
	subsetMasks [][]int
}

func newMultiDP(s []*bio.Sequence) *multiDP {
	return &multiDP{
		seqs:        s,
		table:       nd.NewArray(sizes(s)),
		cached:      nd.NewArray(sizes(s)),
		subsetMasks: generateSubsetMasks(len(s)),
	}
}

/*
generateSubsetMasks takes in the number of sequences being compared
and will generate the column masks which are to be used in the optimal score
computation. these masks are generated once at the start.
ie. generateSubsetMasks(3) will return the following column "masks"
[1 1 1]
[1 1 0]
[1 0 1]
[1 0 0]
[0 1 1]
[0 1 0]
[0 0 1]
notice that [0 0 0] has been removed from the subset masks so as to
prevent alignments from containing columns entirely composed of gaps.
*/
func generateSubsetMasks(numSeqCompared int) [][]int {
	x := []int{}
	for i := 0; i < numSeqCompared; i++ {
		x = append(x, 1)
	}
	return findSubsets(x)
}

// Solve takes in a list of sequences and returns score of the optimal alignment
func Solve(seqStrings []string) float64 {
	mdp := newMultiDP(bio.AsToSeqs(seqStrings))
	return mdp.solve()
}

func (m *multiDP) solve() float64 {
	optScore := m.optimalScore(m.maxIndices())
	// values doubled to be able to use integers during calculations
	// final result is converted to float then divided by two
	// see const block declared above
	return float64(optScore) / 2
}

// uses memoization as opposed to tabulation/dp
// represents optimal score function F(i1, i2, i3, ... , in)
func (m *multiDP) optimalScore(idxs []int) (best int) {
	// base case
	for _, i := range idxs {
		if i <= 0 {
			return
		}
	}
	// have we calculated the score for these indices before?
	if m.cached.At(idxs) == 1 {
		return m.table.At(idxs)
	}
	// see generateSubsetMasks() comment for explaination of subset masks
	// iterate over all possible masks to find the optimal score
	for _, mask := range m.subsetMasks {
		mIdxs, ok := maskedIdxs(idxs, mask)
		// fmt.Printf("mIdxs f(%v): %v - %v\n", idxs, mIdxs, ok)
		if !ok {
			// if an index in the masked indices is below zero, we should skip it
			// because a negative index is invalid and undefined.
			continue
		}
		// find optimal score of masked indices
		optScore := m.optimalScore(mIdxs)

		// maskedBases are the bases (and gaps) given the current indices (idxs)
		// and the mask.
		bases := m.maskedBases(idxs, mask)
		// calculate the score of this column of bases (and gaps) using sum-of-pairs
		score := bio.ColumnSPScore(bases)

		// maintain best score
		best = max(best, optScore+score)
	}
	// save results. mark this specific set of indicies as cached.
	m.table.Set(best, idxs)
	m.cached.Set(1, idxs)
	// fmt.Printf("calced f(%v): %v\n", idxs, best)
	return
}

/////////////////////////
// Helper Functions
/////////////////////////

func sizes(s []*bio.Sequence) (sizes []int) {
	for _, seq := range s {
		sizes = append(sizes, len(seq.Bases)+1)
	}
	return
}

func maskedIdxs(idxs, mask []int) (mIdxs []int, ok bool) {
	for i, idx := range idxs {
		x := idx - mask[i]
		if x < 0 {
			return mIdxs, false
		}
		mIdxs = append(mIdxs, x)
	}
	return mIdxs, true
}

func (m *multiDP) maskedBases(idxs, mask []int) (bases []bio.Base) {
	for i, idx := range idxs {
		var b bio.Base
		if mask[i] == 1 { // not a gap
			b = m.seqs[i].Bases[idx-1]
		} else {
			b = bio.X
		}
		bases = append(bases, b)
	}
	return
}

func (m *multiDP) maxIndices() (indices []int) {
	for _, size := range sizes(m.seqs) {
		indices = append(indices, size-1)
	}
	return
}

func decrementedIndices(idxs []int) (si []int) {
	for _, i := range idxs {
		si = append(si, i-1)
	}
	return
}

func findSubsets(idxs []int) (subsets [][]int) {
	subsets = subsetHelper(idxs, len(idxs)-1)
	return purgeAllGapCase(idxs, subsets)
}

func purgeAllGapCase(idxs []int, subsets [][]int) (ss [][]int) {
	si := decrementedIndices(idxs)
	for _, subset := range subsets {
		if !arraysMatch(subset, si) {
			ss = append(ss, subset)
		}
	}
	return
}

func arraysMatch(s1, s2 []int) bool {
	var count int
	for i, v := range s1 {
		if s2[i] == v {
			count++
		}
	}
	return count == len(s1)
}

func subsetHelper(idxs []int, i int) (subsets [][]int) {
	if i == -1 {
		return append(subsets, []int{})
	}
	x := idxs[i]
	for _, subset := range subsetHelper(idxs, i-1) {
		s1 := cpy(subset)
		s2 := cpy(subset)
		subsets = append(subsets, append(s1, x))
		subsets = append(subsets, append(s2, x-1))
	}
	return
}

func cpy(a []int) (b []int) {
	b = make([]int, len(a))
	copy(b, a)
	return
}

func max(ints ...int) int {
	max := math.MinInt64
	for _, x := range ints {
		if x > max {
			max = x
		}
	}
	return max
}
