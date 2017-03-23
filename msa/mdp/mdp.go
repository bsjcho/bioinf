package mdp

import (
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
)

type multiDP struct {
	seqs  []*Sequence
	table *nd.Array
}

func newMultiDP(s []*Sequence) *multiDP {

	return &multiDP{
		seqs:  s,
		table: nd.NewArray(sizes(s)),
	}
}

func sizes(s []*Sequence) (sizes []int) {
	for _, seq := range s {
		sizes = append(sizes, len(seq.bases))
	}
	return
}

// Solve takes in a list of sequences and returns score of the optimal alignment
func Solve(s []*Sequence) int {
	mdp := newMultiDP(s)
	return mdp.solve()
}

func (m *multiDP) solve() int {
	m.initialize()
	m.fillTables()
	return m.table.At(m.maxIndices())
}

func (m *multiDP) initialize() {
	// TODO
}

func (m *multiDP) fillTables() {
	// TODO
}

func (m *multiDP) maxIndices() (indices []int) {
	for _, size := range sizes(m.seqs) {
		indices = append(indices, size-1)
	}
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
