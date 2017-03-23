package bioinf

import (
	"math"
	"testing"
)

const (
	x1 = "AATTATGG"
)

func TestConversion(t *testing.T) {
	seq := AToSeq(x1)
	rSeq := []Base{A, A, T, T, A, T, G, G}
	for i, b := range seq.Bases {
		if rSeq[i] != b {
			t.Errorf("base conversion failed: %v %v\n", rSeq[i], b)
		}
	}
	math.Abs(2)
}
