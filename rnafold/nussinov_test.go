package rnafold

import "testing"

const (
	seq       = "GGGAAATCC"
	mIR1976   = "GCAGCAAGGAAGGCAGGGGTCCTAAGGTGTGTCCTCCTGCCCTCCTTGCTGT"
	rNA5SP136 = "GTCTGTGGCCATACCACCCAGAACGCACTCGATCTCATCTTATCTCCAAAGCTAAGCATGGTTGGCCTGGTTAGTACTTGGATGGGAGAAAACTGTTATCCTACCTT"
	bCYRN1    = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCTCTCAGGGAGGCTAAGAGGCGGGAGGATAGCTTGAGCCCAGGAGTTCGAGACCTGCCTGGGCAATATAGCGAGACCCCGTTCTCCAGAAAAAGGAAAAAAAAAAACAAAAGACAAAAAAAAAATAAGCGTAACTTCCCTCAAAGCAACAACCCCCCCCCCCCTTT"
)

func TestSolve(t *testing.T) {
	if Solve(seq) != 3 {
		t.Fatal("Solving seq failed.")
	}
	if Solve(mIR1976) != 21 {
		t.Fatal("Solving mIR1976 failed.")
	}
	if Solve(rNA5SP136) != 44 {
		t.Fatal("Solving rNA5SP136 failed.")
	}
	if Solve(bCYRN1) != 73 {
		t.Fatal("Solving bCYRN1 failed.")
	}
}

func BenchmarkMIR1976(b *testing.B) {
	for i := 0; i < b.N; i++ {
		Solve(mIR1976)
	}
}

func BenchmarkRNA5SP136(b *testing.B) {
	for i := 0; i < b.N; i++ {
		Solve(rNA5SP136)
	}
}

func BenchmarkBCYRN1(b *testing.B) {
	for i := 0; i < b.N; i++ {
		Solve(bCYRN1)
	}
}
