package rnafold

import (
	bio "github.com/bsjcho/bioinf"
)

// RNA secondary structure using the Nussinov algorithm.
// Input:
// - A sequence composed of nucleotide bases (A,C,G,T)
// Output:
// - The number of base pair matchings

type nussinov struct {
	sequence string
	matrix   [][]int
}

// FoldScore performs the nussinov algorithm on a sequence, returning the max pairs.
func FoldScore(seq string) (pairCount int) {
	n := &nussinov{sequence: seq}
	initialize(n)
	return foldScore(0, len(n.sequence)-1, n)
}

func initialize(n *nussinov) {
	seqLen := len(n.sequence)
	n.matrix = bio.Slice2D(seqLen, seqLen, -1)
}

func foldScore(i, j int, n *nussinov) int {
	if j-i < 2 {
		return 0
	}
	if n.matrix[i][j] == -1 {
		n.matrix[i][j] = bio.Max(
			foldScore(i+1, j, n),
			foldScore(i, j-1, n),
			foldScore(i+1, j-1, n)+matchScore(i, j, n),
		)
		for k := i + 1; k < j; k++ {
			n.matrix[i][j] = bio.Max(
				n.matrix[i][j],
				foldScore(i, k, n)+foldScore(k+1, j, n))
		}
	}
	return n.matrix[i][j]
}

func matchScore(i, j int, n *nussinov) int {
	ib := string(n.sequence[i])
	jb := string(n.sequence[j])
	t := ib + jb
	res := t == "AU" || t == "CG" || t == "GU" ||
		t == "UA" || t == "GC" || t == "UG"
	if res {
		return 1
	}
	return 0
}
