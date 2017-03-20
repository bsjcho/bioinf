package nussinov

import "testing"

const (
	mIR1976   = "GCAGCAAGGAAGGCAGGGGTCCTAAGGTGTGTCCTCCTGCCCTCCTTGCTGT"
	rNA5SP136 = "GTCTGTGGCCATACCACCCAGAACGCACTCGATCTCATCTTATCTCCAAAGCTAAGCATGGTTGGCCTGGTTAGTACTTGGATGGGAGAAAACTGTTATCCTACCTT"
	bCYRN1    = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCTCTCAGGGAGGCTAAGAGGCGGGAGGATAGCTTGAGCCCAGGAGTTCGAGACCTGCCTGGGCAATATAGCGAGACCCCGTTCTCCAGAAAAAGGAAAAAAAAAAACAAAAGACAAAAAAAAAATAAGCGTAACTTCCCTCAAAGCAACAACCCCCCCCCCCCTTT"
)

func TestSolve(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping test in short mode.")
	}
	seq := "GGGAAAUCC"
	t.Logf("When considering the sequence %v\n", seq)
	pairCount := Solve(seq)
	if pairCount == 3 {
		t.Logf("\tShould be able to calculate max base pairs\n")
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
