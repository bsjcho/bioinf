package nussinov

import "testing"

const (
	checkMark = "\u2713"
	ballotX   = "\u2717"
)

func TestSolve(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping test in short mode.")
	}
	seq := "GGGAAAUCC"
	t.Logf("When consiering the sequence %v\n", seq)
	pairCount := Solve(seq)
	if pairCount == 3 {
		t.Logf("\tShould be able to calculate max base pairs %v\n", checkMark)
	}
}
