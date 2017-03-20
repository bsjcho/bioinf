package rnafold

import "math"

// RNA secondary structure using the Nussinov algorithm.
// Input:
// - A sequence composed of nucleotide bases (A,C,G,T)
// Output:
// - The number of base pair matchings

type nussinov struct {
	sequence string
	matrix   [][]int
}

// Solve performs the nussinov algorithm on a sequence, returning the max pairs.
func Solve(seq string) (pairCount int) {
	n := &nussinov{sequence: seq}
	initialize(n)
	return maxPairs(0, len(n.sequence)-1, n)
}

func initialize(n *nussinov) {
	seqLen := len(n.sequence)
	n.matrix = slice2D(seqLen, seqLen, -1)
}

func maxPairs(i, j int, n *nussinov) int {
	if i == j || j == i-1 {
		return 0
	}
	if n.matrix[i][j] > -1 {
		return n.matrix[i][j]
	}
	n.matrix[i][j] = max(
		maxPairs(i+1, j, n),
		maxPairs(i, j-1, n),
		maxPairs(i+1, j-1, n)+score(i, j, n),
	)
	for k := i + 1; k < j; k++ {
		n.matrix[i][j] = max(
			n.matrix[i][j],
			maxPairs(i, k, n)+maxPairs(k+1, j, n))
	}
	return n.matrix[i][j]
}

func score(i, j int, n *nussinov) int {
	ib := string(n.sequence[i])
	jb := string(n.sequence[j])
	if ib == "G" && jb == "C" ||
		ib == "C" && jb == "G" ||
		ib == "A" && (jb == "T" || jb == "U") ||
		(ib == "T" || ib == "U") && jb == "A" {
		return 1
	}
	return 0
}

func slice2D(rows, cols, defVal int) [][]int {
	m := make([][]int, rows)
	x := make([]int, rows*cols)
	for j := range x {
		x[j] = defVal
	}
	for i := range m {
		m[i], x = x[:cols], x[cols:]
	}
	return m
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
