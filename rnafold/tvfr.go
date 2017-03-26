package rnafold

import bio "github.com/bsjcho/bioinf"

type tvfrNussinov struct {
	n        int
	sequence string
	matrix   [][]int
	score    int
	qSq      int
	q        int
	R        []int
}

// TVFRFoldScore returns the fold score using nussinov
// with the two vector four russians speedup.
func TVFRFoldScore(seq string) int {
	var tvfr tvfrNussinov
	tvfr.initialize(seq)
	// TODO
	return 0
}

func (t *tvfrNussinov) initialize(seq string) {
	t.n = len(seq)
	t.sequence = seq
	t.matrix = bio.Slice2D(t.n, t.n, 0)
	// groupSize := math.Log(x) log10(x)/log10(5)
}
