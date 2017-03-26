package rnafold

import "testing"

const (
	mIR1976   = "GCAGCAAGGAAGGCAGGGGUCCUAAGGUGUGUCCUCCUGCCCUCCUUGCUGU"
	rNA5SP136 = "GUCUGUGGCCAUACCACCCAGAACGCACUCGAUCUCAUCUUAUCUCCAAAGCUAAGCAUGGUUGGCCUGGUUAGUACUUGGAUGGGAGAAAACUGUUAUCCUACCUU"
	bCYRN1    = "GGCCGGGCGCGGUGGCUCACGCCUGUAAUCCCAGCUCUCAGGGAGGCUAAGAGGCGGGAGGAUAGCUUGAGCCCAGGAGUUCGAGACCUGCCUGGGCAAUAUAGCGAGACCCCGUUCUCCAGAAAAAGGAAAAAAAAAAACAAAAGACAAAAAAAAAAUAAGCGUAACUUCCCUCAAAGCAACAACCCCCCCCCCCCUUU"
)

func TestSolve(t *testing.T) {
	t1 := Solve(mIR1976)
	t.Log(t1)
	if t1 != 23 {
		t.Fatal("Solving mIR1976 failed.")
	}
	t2 := Solve(rNA5SP136)
	t.Log(t2)
	if t2 != 42 {
		t.Fatal("Solving rNA5SP136 failed.")
	}
	t3 := Solve(bCYRN1)
	t.Log(t3)
	if t3 != 69 {
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
