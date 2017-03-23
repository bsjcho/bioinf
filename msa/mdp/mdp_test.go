package mdp

import (
	"fmt"
	"testing"
)

const (
	x1 = "AATTATGG"
	x2 = "ACATTGTTG"
	x3 = "GCCAGGAGG"
	x4 = "AATTTTGAGG"

	x5 = "AA"
	x6 = "AA"
	x7 = "AA"
	x8 = "AA"
)

func TestFindSubsets(t *testing.T) {
	idxs := []int{1, 1, 1}
	ss := findSubsets(idxs)
	fmt.Println(len(ss))
	for _, z := range ss {
		fmt.Println(z)
	}
}

func TestGenSubMasks(t *testing.T) {
	ss := generateSubsetMasks(4)
	for _, z := range ss {
		fmt.Println(z)
	}
}

func TestConversion(t *testing.T) {
	seq := convertStringSequence(x1)
	rSeq := []Base{A, A, T, T, A, T, G, G}
	for i, b := range seq.bases {
		if rSeq[i] != b {
			t.Errorf("base conversion failed: %v %v\n", rSeq[i], b)
		}
	}
}

func solveSeqs(seqStrings []string) float64 {
	seqs := []*Sequence{}
	for _, seqStr := range seqStrings {
		seqs = append(seqs, convertStringSequence(seqStr))
	}
	mdp := newMultiDP(seqs)
	return mdp.solve()
}

func TestMDP(t *testing.T) {
	fmt.Println(solveSeqs([]string{x1, x2, x3, x4}))
}

func TestMDPSimple(t *testing.T) {
	optScore := solveSeqs([]string{x5, x6, x7, x8})
	fmt.Println(optScore)
	if optScore != 36 {
		t.Error("wrong score for simple MDP")
	}
}
