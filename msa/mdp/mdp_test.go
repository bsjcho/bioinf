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

func TestMDP(t *testing.T) {
	optScore := Solve([]string{x1, x2, x3, x4})
	t.Log(optScore)
	if optScore != 45 {
		t.Error("Incorrect score.")
	}
}

func TestMDPSimple(t *testing.T) {
	optScore := Solve([]string{x5, x6, x7, x8})
	t.Log(optScore)
	if optScore != 36 {
		t.Error("Incorrect score.")
	}
}
