package mdp

import (
	"fmt"
	"math"

	"github.com/bsjcho/nd"
)

// Sequence represents a nucleotide base sequence
type Sequence struct {
	bases []Base
}

// NewSequence is a Sequence constructor
func NewSequence() *Sequence {
	return &Sequence{bases: []Base{}}
}

// Base represents a nucleotide base
type Base int

// A ... enum represents a nucleotide
const (
	A Base = iota
	C
	G
	T
	X // represents a gap "-"
)

func convertStringSequence(seq string) *Sequence {
	s := NewSequence()
	for _, b := range seq {
		s.bases = append(s.bases, convertStringBase(string(b)))
	}
	return s
}

func convertStringBase(b string) Base {
	switch b {
	case "A":
		return A
	case "C":
		return C
	case "G":
		return G
	case "T":
		return T
	default:
		return X
	}
}

type multiDP struct {
	seqs        []*Sequence
	table       *nd.Array
	cached      *nd.Array
	subsetMasks [][]int
}

func newMultiDP(s []*Sequence) *multiDP {

	return &multiDP{
		seqs:        s,
		table:       nd.NewArray(sizes(s)),
		cached:      nd.NewArray(sizes(s)),
		subsetMasks: generateSubsetMasks(len(s)),
	}
}

func generateSubsetMasks(size int) [][]int {
	x := []int{}
	for i := 0; i < size; i++ {
		x = append(x, 1)
	}
	return findSubsets(x)
}

func sizes(s []*Sequence) (sizes []int) {
	for _, seq := range s {
		sizes = append(sizes, len(seq.bases)+1)
	}
	return
}

// Solve takes in a list of sequences and returns score of the optimal alignment
func Solve(s []*Sequence) float64 {
	mdp := newMultiDP(s)
	return mdp.solve()
}

func (m *multiDP) solve() float64 {
	maxIdxs := m.maxIndices()
	fmt.Println(maxIdxs)
	return float64(m.f(maxIdxs)) / 2
}

func (m *multiDP) f(idxs []int) (best int) {
	for _, i := range idxs {
		if i <= 0 {
			return
		}
	}
	if m.cached.At(idxs) == 1 {
		return m.table.At(idxs)
	}
	for _, mask := range m.subsetMasks {
		mIdxs, ok := maskedIdxs(idxs, mask)
		// fmt.Printf("mIdxs f(%v): %v - %v\n", idxs, mIdxs, ok)
		if !ok { // if an idx is below 0, don't consider this option
			continue
		}
		bestf := m.f(mIdxs)
		score := m.score(m.maskedBases(idxs, mask))
		best = max(best, bestf+score)
	}
	m.table.Set(best, idxs)
	m.cached.Set(1, idxs)
	// fmt.Printf("calced f(%v): %v\n", idxs, best)
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

func (m *multiDP) maskedBases(idxs, mask []int) (bases []Base) {
	for i, idx := range idxs {
		var b Base
		if mask[i] == 1 { // not a gap
			b = m.seqs[i].bases[idx-1]
		} else {
			b = X
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

func (m *multiDP) score(bases []Base) (sum int) {
	for i, bi := range bases[:len(bases)-1] {
		for _, bj := range bases[i+1:] {
			sum += pairScore(bi, bj)
		}
	}
	return
}

const (
	// values doubled to be able to use integers during calculations
	match    = 6
	mismatch = -4
	gap      = -3
)

func pairScore(b1, b2 Base) int {
	if b1 == X && b2 == X {
		return 0
	}
	if (b1 == X && b2 != X) ||
		(b2 == X && b1 != X) {
		return gap
	}
	if b1 != b2 {
		return mismatch
	}
	return match
}

///////////////////////// Helpers

func subbedIndices(idxs []int) (si []int) {
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
	si := subbedIndices(idxs)
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
